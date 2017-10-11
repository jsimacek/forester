/*
 * Copyright (C) 2011 Jiri Simacek
 *
 * This file is part of forester.
 *
 * forester is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * forester is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with forester.  If not, see <http://www.gnu.org/licenses/>.
 */

// Standard library headers
#include <ostream>

// Code Listener headers
#include <cl/storage.hh>

// Forester headers
#include "abstraction.hh"
#include "config.h"
#include "executionmanager.hh"
#include "fixpoint.hh"
#include "folding.hh"
#include "forestautext.hh"
#include "normalization.hh"
#include "regdef.hh"
#include "splitting.hh"
#include "streams.hh"
#include "ufae.hh"
#include "utils.hh"
#include "virtualmachine.hh"

// anonymous namespace
namespace
{
struct ExactTMatchF
{
	bool operator()(
		const TreeAut::Transition&  t1,
		const TreeAut::Transition&  t2)
	{
		return TreeAut::GetSymbol(t1) == TreeAut::GetSymbol(t2);
	}
};

std::vector<size_t> findSelectors(
		const FAE&                    fae,
		const TreeAut::Transition&    t)
{
	std::vector<size_t> index;

    for (const AbstractBox* aBox : TreeAut::GetSymbol(t)->getNode())
    {
        if (aBox->isBox())
        {
            index.push_back(aBox->getOrder());
        }
	}

	return index;
}

struct SmartTMatchF
{
public:
	SmartTMatchF()
	{
	}

	bool operator()(
		const TreeAut::Transition&  t1,
		const TreeAut::Transition&  t2)
	{
        label_type l1 = TreeAut::GetSymbol(t1);
        label_type l2 = TreeAut::GetSymbol(t2);
		if (l1->isNode() && l2->isNode())
		{
			if (!FA_ALLOW_STACK_FRAME_ABSTRACTION)
			{
				if ((static_cast<const TypeBox *>(l1->nodeLookup(-1).aBox))->getName().find("__@") == 0)
				{
					return l1 == l2;
				}
			}

			return l1->getTag() == l2->getTag();
		}

		return l1 == l2;
	}
};

class SmarterTMatchF
{
private:  // data members

	const FAE& fae_;

public:   // methods

	SmarterTMatchF(const FAE& fae) :
		fae_(fae)
	{ }

	bool operator()(
		const TreeAut::Transition&  t1,
		const TreeAut::Transition&  t2)
	{
        label_type l1 = TreeAut::GetSymbol(t1);
        label_type l2 = TreeAut::GetSymbol(t2);
		if (!l1->isNode() || !l2->isNode())
			return l1 == l2;

		if (l1->getTag() != l2->getTag())
			return false;

		if (&t1.GetChildren() == &t2.GetChildren())
			return true;

		if (t1.GetChildrenSize() != t2.GetChildrenSize())
			return false;

		for (size_t i = 0; i < t1.GetChildrenSize(); ++i)
		{
			size_t s1 = t1.GetNthChildren(i), s2 = t2.GetNthChildren(i), ref;

			if (s1 == s2)
				continue;

			if (FA::isData(s1))
			{
				if (!this->fae_.getRef(s1, ref))
					return false;

				if (FA::isData(s2))
					return false;
			} else
			{
				if (FA::isData(s2) && !fae_.getRef(s2, ref))
					return false;
			}
		}

		return true;
	}
};

struct CompareVariablesF
{
	bool operator()(
		const FAE&            fae,
		size_t                i,
		const TreeAut&        ta1,
		const TreeAut&        ta2)
	{
		VirtualMachine vm(fae);

		bool isFixedComp = true;
		for (size_t j = 0; j < FIXED_REG_COUNT; ++j)
		{
			if (i == vm.varGet(j).d_ref.root)
			{
				isFixedComp = false;
				break;
			}
		}

		if (isFixedComp)
		{
			return true;
		}

		const TreeAut::Transition& t1 = ta1.getAcceptingTransition();
		const TreeAut::Transition& t2 = ta2.getAcceptingTransition();
		return (TreeAut::GetSymbol(t1) == TreeAut::GetSymbol(t2)) &&
                (t1.GetChildren() == t2.GetChildren());
	}
};

struct FuseNonFixedF
{
	bool operator()(size_t root, const FAE* fae)
	{
		// Preconditions
		assert(nullptr != fae);

		VirtualMachine vm(*fae);

		for (size_t i = 0; i < FIXED_REG_COUNT; ++i)
		{
			if (root == vm.varGet(i).d_ref.root)
			{
				return false;
			}
		}

		return true;
	}
};


/**
 * @brief  Reorders components into the canonical order
 *
 * This function reorders the components of the FA @p fae @e without @e merging
 * them together.
 *
 * @param[in]      state  The state of the symbolic execution
 * @param[in,out]  fae    The forest automaton to be reordered
 */
void reorder(
	const SymState*   state,
	FAE&              fae)
{
	fae.unreachableFree();
	Normalization::normalizeWithoutMerging(fae, state);

	FA_DEBUG_AT(3, "after reordering: " << std::endl << fae);
}

bool testInclusion(
	FAE&                           fae,
	TreeAut&                       fwdConf,
	UFAE&                          fwdConfWrapper)
{
	TreeAut ta = TreeAut::createTAWithSameTransitions(fwdConf);

	Index<size_t> index;

	fae.unreachableFree();

	fwdConfWrapper.fae2ta(ta, index, fae);

	if (TreeAut::subseteq(ta, fwdConf))
		return true;

	fwdConfWrapper.join(ta, index);

	ta.clear();

	fwdConf.minimized(ta);
	fwdConf = ta;

	return false;
}

struct CopyNonZeroRhsF
{
	bool operator()(const TreeAut::Transition* transition) const
	{
		return transition->GetParent() != 0;
	}
};


void getCandidates(
	std::set<size_t>&               candidates,
	const FAE&                      fae)
{
	std::unordered_map<
		std::vector<std::pair<int, size_t>>,
		std::set<size_t>,
		boost::hash<std::vector<std::pair<int, size_t>>>
	> partition;

	for (size_t i = 0; i < fae.getRootCount(); ++i)
	{
		std::vector<std::pair<int, size_t>> tmp;

		fae.connectionGraph.getRelativeSignature(tmp, i);

		partition.insert(std::make_pair(tmp, std::set<size_t>())).first->second.insert(i);
	}

	candidates.clear();

	for (auto& tmp : partition)
	{
		if (tmp.second.size() > 1)
			candidates.insert(tmp.second.begin(), tmp.second.end());
	}
}
} // namespace

void FixpointBase::initFoldedRoots()
{
	abstrIteration_ = 0;
	iterationToFoldedRoots_.clear();
}


void FI_abs::preFold(Folding::StateToBoxInfoMap& stateToBoxInfoMap, const std::shared_ptr<FAE>& src)
{
	FAE fae(*src);

	FA_DEBUG_AT(3, "before pre-folding: " << std::endl << fae);

	std::set<size_t> forbidden = Normalization::computeForbiddenSet(fae);

	size_t level = 2;

#if FA_BOX_APPROXIMATION
	auto foldResult = Folding::fold(fae, boxMan_, boxAntichain_, forbidden, level, true, true, stateToBoxInfoMap).second;
#else
	auto foldResult = Folding::fold(fae, boxMan_, forbidden, level, true, true, stateToBoxInfoMap).second;
#endif

	FAE old(fae, boxMan_);

	do
	{
		forbidden = Normalization::computeForbiddenSet(fae);

		Normalization::normalize(fae, nullptr, forbidden, false);

		abstract(fae);

		forbidden.clear();
		for (size_t i = 0; i < FIXED_REG_COUNT; ++i)
			forbidden.insert(VirtualMachine(fae).varGet(i).d_ref.root);

		old = fae;

#if FA_BOX_APPROXIMATION
		foldResult = Folding::fold(fae, boxMan_, boxAntichain_, forbidden, level++, false, true, stateToBoxInfoMap).second;
#else
		foldResult = Folding::fold(fae, boxMan_, forbidden, level++, false, true, stateToBoxInfoMap).second;
#endif
	}
	while (foldResult && !FAE::subseteq(fae, old));

	FA_DEBUG_AT(3, "after pre-folding: " << std::endl << fae);

	for (size_t root = 0; root < fae.getRootCount(); ++root)
		Folding::analyzeBoxes(stateToBoxInfoMap, *fae.getRoot(root));

	for (auto& stateBoxInfoPair : stateToBoxInfoMap)
	{
		FA_DEBUG_AT(3, "box analysis at state: " << stateBoxInfoPair.first);
		for (auto& boxBoolPair : stateBoxInfoPair.second)
			FA_DEBUG_AT(3, boxBoolPair.first->getName() << ": " << boxBoolPair.second);
	}
}


bool FixpointBase::fold(
	const std::shared_ptr<FAE>&       fae,
	std::set<size_t>&                 forbidden,
	size_t				  level,
	bool				  discover3only,
	Folding::StateToBoxInfoMap&	  stateToBoxInfoMap)
{
#if FA_BOX_APPROXIMATION
	auto foldResult = Folding::fold(*fae, boxMan_, boxAntichain_, forbidden, level, discover3only, false, stateToBoxInfoMap);
#else
	auto foldResult = Folding::fold(*fae, boxMan_, forbidden, level, discover3only, false, stateToBoxInfoMap);
#endif
	iterationToFoldedRoots_[abstrIteration_] = foldResult.first;
	++abstrIteration_;

	return foldResult.second;
}


SymState* FixpointBase::reverseAndIsect(
	ExecutionManager&                      execMan,
	const SymState&                        fwdPred,
	const SymState&                        bwdSucc) const
{
	(void)fwdPred;
	SymState* tmpState = execMan.copyState(bwdSucc);

	assert(iterationToFoldedRoots_.size() - faeAtIteration_.size() <= 1);
	for (int i = abstrIteration_ - 1; i >= 0; --i)
	{
		assert(i < iterationToFoldedRoots_.size());

		for (const auto& rootToBoxes : iterationToFoldedRoots_.at(i))
		{
			for (const auto& boxes : rootToBoxes.second)
			{
				std::shared_ptr<FAE> fae = std::shared_ptr<FAE>(new FAE(*(tmpState->GetFAE())));
				if (fae->getRootCount() == 0)
				{
					return tmpState;
				}

				std::vector<FAE*> unfolded;
				Splitting splitting(*fae);
				splitting.isolateSet(
						unfolded,
						rootToBoxes.first,
						0,
						findSelectors(*fae, fae->getRoot(rootToBoxes.first)->getAcceptingTransition()));
						//fae->getType(rootToBoxes.first)->getSelectors());
				assert(unfolded.size() == 1);
				fae = std::shared_ptr<FAE>(unfolded.at(0));
				try
				{
					Unfolding(*fae).unfoldBox(rootToBoxes.first, boxes);
				}
				catch (std::runtime_error e)
				{
					;
				}
				tmpState->SetFAE(fae);
			}
		}

		if (faeAtIteration_.count(i))
		{
			SymState st;
			st.init(*tmpState);
			st.SetFAE(faeAtIteration_.at(i));
			tmpState->Intersect(st);
		}
	}

	// perform intersection
	tmpState->Intersect(fwdPred);

	FA_DEBUG_AT(1, "Executing !!VERY!! suspicious reverse operation FixpointBase");
	return tmpState;
}

void FixpointBase::fixpoint(ExecutionManager& execMan, SymState& state, std::shared_ptr<FAE> fae)
{
	fae->updateConnectionGraph();

	Normalization::normalize(*fae, &state, Normalization::computeForbiddenSet(*fae), true);

	if (testInclusion(*fae, fwdConf_, fwdConfWrapper_))
	{
		FA_DEBUG_AT(3, "hit");

		execMan.pathFinished(&state);
	} else
	{
		FA_DEBUG_AT(1, "extending fixpoint\n" << *fae);

		SymState* tmpState = execMan.createChildState(state, next_);
		tmpState->SetFAE(fae);

		execMan.enqueue(tmpState);
	}
}

void FI_abs::abstract(
	FAE&                 fae)
{
	fae.unreachableFree();

	FA_DEBUG_AT(3, "before abstraction: " << std::endl << fae);

	if (FA_FUSION_ENABLED)
	{
		// merge fixpoint
		std::vector<FAE*> tmp;

		ContainerGuard<std::vector<FAE*>> g(tmp);

		FAE::loadCompatibleFAs(
			/* the result */ tmp,
			fwdConf_,
			ta_,
			boxMan_,
			fae,
			0,
			CompareVariablesF()
		);

		for (size_t i = 0; i < tmp.size(); ++i)
		{
			FA_DEBUG_AT(3, "accelerator " << std::endl << *tmp[i]);
		}

		fae.fuse(tmp, FuseNonFixedF());
		FA_DEBUG_AT(3, "fused " << std::endl << fae);
	}

	// abstract
	Abstraction abstraction(fae);

	if ((FA_START_WITH_PREDICATE_ABSTRACTION || !predicates_.empty()) && FA_USE_PREDICATE_ABSTRACTION)
	{	// for predicate abstraction
		abstraction.predicateAbstraction(this->getPredicates());
	}
	else
	{	// for finite height abstraction
		FA_DEBUG_AT(1,"Height abstraction of height " << FA_ABS_HEIGHT);

		// the roots that will be excluded from abstraction
		std::vector<bool> excludedRoots(fae.getRootCount(), false);
		for (size_t i = 0; i < FIXED_REG_COUNT; ++i)
		{
			excludedRoots[VirtualMachine(fae).varGet(i).d_ref.root] = true;
		}

		for (size_t i = 0; i < fae.getRootCount(); ++i)
		{
			if (!excludedRoots[i])
			{
				abstraction.heightAbstraction(i, FA_ABS_HEIGHT, SmartTMatchF());

//				abstraction.heightAbstraction(i, FA_ABS_HEIGHT, SmarterTMatchF(fae));
			}
		}
	}

	FA_DEBUG_AT(3, "after abstraction: " << std::endl << fae);
}


// FI_abs
void FI_abs::execute(ExecutionManager& execMan, SymState& state)
{
	this->initFoldedRoots();

	std::shared_ptr<FAE> fae = std::shared_ptr<FAE>(new FAE(*(state.GetFAE())));

	fae->updateConnectionGraph();

	std::set<size_t> forbidden = Normalization::computeForbiddenSet(*fae);

	Normalization::normalizeWithoutMerging(*fae, &state, forbidden, false);

	Folding::StateToBoxInfoMap stateToBoxInfoMap;

#if FA_ALLOW_FOLDING
//	this->preFold(stateToBoxInfoMap, fae);

	this->fold(fae, forbidden, 2, true, stateToBoxInfoMap);
#endif

	FAE old(*fae, boxMan_);

	size_t level = 2;

	do
	{
		forbidden = Normalization::computeForbiddenSet(*fae);

		Normalization::normalize(*fae, &state, forbidden, false);

		faeAtIteration_[abstrIteration_] = std::shared_ptr<FAE>(new FAE(*fae));

		abstract(*fae);

		forbidden.clear();
		for (size_t i = 0; i < FIXED_REG_COUNT; ++i)
			forbidden.insert(VirtualMachine(*fae).varGet(i).d_ref.root);

		old = *fae;
	}
#if FA_ALLOW_FOLDING
	while (this->fold(fae, forbidden, level++, false, stateToBoxInfoMap) && !FAE::subseteq(*fae, old));
#else
	while (false);
#endif

	FA_DEBUG_AT(3, "after folding: " << std::endl << *fae);

	for (size_t root = 0; root < fae->getRootCount(); ++root)
	{
		if (fae->getRoot(root) != nullptr)
			Unfolding(*fae).unfoldSingletons(root);
	}

	FA_DEBUG_AT(3, "after unfolding: " << std::endl << *fae);

//	FA_DEBUG_AT_MSG(1, &this->insn()->loc, "AbsInt end " << *fae);

	fixpoint(execMan, state, fae);
}

// FI_fix
void FI_fix::execute(ExecutionManager& execMan, SymState& state)
{
	fixpoint(execMan, state, std::shared_ptr<FAE>(new FAE(*(state.GetFAE()))));
}

void FI_abs::printDebugInfoAboutPredicates() const
{
    FA_DEBUG_AT(1, "Number of predicates " <<
            this->getPredicates().size() << " of " << this->insn());
    std::ostringstream os;
    this->printPredicates(os);
    FA_DEBUG_AT(1, os.str());
}

void FI_abs::printPredicates(std::ostringstream& os) const
{
	for (const std::shared_ptr<const TreeAut>& pred : this->getPredicates())
	{
		os << "\n---------------------------------------------------\n"
			<< *this << this->insn();

		const CodeStorage::Insn* clInsn = this->insn();
		assert(nullptr != clInsn);
		assert(nullptr != clInsn->bb);
		if (clInsn->bb->front() == clInsn)
		{
			os << " (" << clInsn->bb->name() << ")";
		}

		if (pred != nullptr)
		{
			os << ": " << *this << "\n" << *pred << '\n';
		}
		else
		{
			os << ": " << *this << "\n" << "nullptr" << '\n';
		}
	}
}
