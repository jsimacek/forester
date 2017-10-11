/*
 * Copyright (C) 2013 Jiri Simacek
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

// Forester headers
#include "folding.hh"

#include "config.h"
#include "normalization.hh"

namespace
{	// anonymous namespace

/**
 * @brief  Adds a mapping of a cutpoint to selector into a map
 *
 * This function updates a selector map, which maps a cutpoint number to the
 * offset of the lowest selector that leads to the given cutpoint, with a new
 * pair of a cutpoint and a selector offset (or updates the image of an
 * already existing cutpoint).
 *
 * @param[in,out]  m          The cutpoint to selector offset map that is to
 *                            be updated
 * @param[in]      selector   The selector offset
 * @param[in]      signature  Signature of the output of the selector
 */
void updateSelectorMap(
	std::unordered_map<size_t, size_t>&          m,
	size_t                                       selector,
	const ConnectionGraph::CutpointSignature&    signature)
{
	for (const ConnectionGraph::CutpointInfo& cutpoint : signature)
	{	// for all cutpoints in the signature

		// insert the mapping 'cutpoint' -> 'selector offset'
		auto p = m.insert(std::make_pair(cutpoint.root, selector));

		if (!p.second && (p.first->second > selector))
		{	// in case the mapping is there and the original 'selector offset' is
			// bigger than the new one, set the new (smaller) one as the image of
			// 'cutpoint'
			p.first->second = selector;
		}
	}
}


/**
 * @brief  Computes cutpoint-to-selector mapping
 *
 * This function computes for a given transition @p t the cutpoint-to-selector
 * mapping, i.e. a map telling for a given cutpoint which selector (or, more
 * precisely, the selector at which @e offset) one needs to take in order to
 * reach the cutpoint.
 *
 * @param[out]  selectorMap  The resulting cutpoint-to-selector map
 * @param[in]   t            The transition
 * @param[in]   stateMap     The mapping of states to their signatures
 */
void computeSelectorMapTrans(
	std::unordered_map<size_t, size_t>&                    selectorMap,
	const TreeAut::Transition&                             t,
	const ConnectionGraph::StateToCutpointSignatureMap&    stateMap)
{
	size_t lhsOffset = 0;

	for (const AbstractBox* absBox : TreeAut::GetSymbol(t)->getNode())
	{	// for all boxes in the transition
		switch (absBox->getType())
		{
			case box_type_e::bSel:
			{	// for ordinary selectors 

				// find the cutpoint signature of the state at position 'lhsOffset' in
				// the transition
				auto iter = stateMap.find(t.GetNthChildren(lhsOffset));

				assert(iter != stateMap.end());

				updateSelectorMap(
					selectorMap,
					(static_cast<const SelBox*>(absBox))->getData().offset,
					iter->second
				);

				break;
			}

			case box_type_e::bBox:
			{	// for nested FA
				const Box* box = static_cast<const Box*>(absBox);

				for (size_t i = 0; i < box->getArity(); ++i)
				{
					auto iter = stateMap.find(t.GetNthChildren(lhsOffset + i));

					assert(iter != stateMap.end());

					updateSelectorMap(
						selectorMap, box->getSelector(i), iter->second
					);
				}

				break;
			}

			default: break;
		}

		lhsOffset += absBox->getArity();
	}
}


/**
 * @brief  Checks whether two cutpoint signatures are compatible
 *
 * Checks whether two cutpoint signatures are compatible.
 *
 * @param[in]  s1  The LHS cutpoint signature
 * @param[in]  s2  The RHS cutpoint signature
 *
 * @returns  @p true of @p s1 and @p s2 are compatible, @p false otherwise
 */
bool isSignaturesCompatible(
	const ConnectionGraph::CutpointSignature&    s1,
	const ConnectionGraph::CutpointSignature&    s2)
{
	if (s1.size() != s2.size())
		return false;

	for (size_t i = 0; i < s1.size(); ++i)
	{
		if (s1[i].root != s2[i].root)
			return false;
/*
		if (*s1[i].fwdSelectors.begin() != *s2[i].fwdSelectors.begin())
			return false;
*/
		if (s1[i].bwdSelector != s2[i].bwdSelector)
			return false;

		if (s1[i].defines != s2[i].defines)
			return false;
	}

	return true;
}


/**
 * @brief  Extracts a selector leading to given cutpoint
 *
 * This method checks a cutpoint-to-selector map @p selectorMap to find the
 * (lowest offset) selector leading to the cutpoint @p target
 *
 * @param[in]  selectorMap  The cutpoint-to-selector map
 * @param[in]  target       Index of the target cutpoint
 *
 * @returns  Offset of the lowest selector leading to cutpoint @p target
 */
size_t extractSelector(
	const std::unordered_map<size_t, size_t>&    selectorMap,
	size_t                                       target)
{
	auto iter = selectorMap.find(target);

	assert(iter != selectorMap.end());

	return iter->second;
}


/**
 * @brief  Creates the mapping of indices of components in the box to selectors
 *
 * This function creates the mapping of indices of components in the box to
 * selectors, i.e. it maps a component number to offset of the selector that
 * maps to the component.
 *
 * @param[in]  selectorMap  The cutpoint-to-selector map
 * @param[in]  root         Index of the input tree automaton
 * @param[in]  index        The map of each cutpoint to the order in which it is
 *                          referenced in the box (or to -1 if it is not referenced)
 *
 * @returns  The map of indices of components to offsets of corresponding selectors
 */
std::vector<size_t> extractInputMap(
	const std::unordered_map<size_t, size_t>&    selectorMap,
	size_t                                       root,
	const std::vector<size_t>&                   index)
{
	// Preconditions
	assert(0 == index[root]);

	std::vector<size_t> inputMap(selectorMap.size());

	size_t count = 0;

	for (const std::pair<size_t, size_t>& cutpointSelectorPair : selectorMap)
	{	// for all cutpoint-to-selector pairs
		if (cutpointSelectorPair.first == root)
		{ // reference to root does not appear in the box interface
			continue;
		}

		// make sure that the cutpoint has a valid number
		assert(cutpointSelectorPair.first < index.size());

		if (index[cutpointSelectorPair.first] == static_cast<size_t>(-1))
		{	// in the case the cutpoint is not referenced in the box
			continue;
		}

		assert(index[cutpointSelectorPair.first] >= 1);
		assert(index[cutpointSelectorPair.first] < inputMap.size() + 1);

		// we now say that the output component at given index is mapped to the
		// given selector
		inputMap[index[cutpointSelectorPair.first] - 1] = cutpointSelectorPair.second;

		++count;
	}

	inputMap.resize(count);

	return inputMap;
}

} // namespace


std::pair<Folding::TreeAutShPtr, Folding::TreeAutShPtr> Folding::separateCutpoint(
	ConnectionGraph::CutpointSignature&            boxSignature,
	const size_t                                   root,
	const size_t                                   state,
	const size_t                                   cutpoint)
{
	TreeAutShPtr ta  = TreeAutShPtr(fae_.allocTA());
	TreeAutShPtr tmp = TreeAutShPtr(fae_.allocTA());

	this->componentCut(*ta, *tmp, boxSignature, root, state, cutpoint);

	TreeAutShPtr tmp2 = TreeAutShPtr(fae_.allocTA());

	tmp->unreachableFree(*tmp2);

	return std::make_pair(ta, tmp2);
}


const ConnectionGraph::StateToCutpointSignatureMap& Folding::getSignatures(
	size_t      root)
{
	// Preconditions
	assert(root < signatureMap_.size());

	if (!signatureMap_[root].first)
	{	// if the signature is not valid, recompute it
		ConnectionGraph::computeSignatures(
			signatureMap_[root].second, *fae_.getRoot(root)
		);

		signatureMap_[root].first = true;
	}

	assert(signatureMap_[root].first);

	return signatureMap_[root].second;
}


bool Folding::discover1(
	size_t                       root,
	const std::set<size_t>&      forbidden,
	size_t			     level,
	bool                         analysis,
	StateToBoxInfoMap& 	     stateToBoxInfoMap,
	std::vector<const Box *>*     discoveredBoxes)
{
	// Preconditions
	assert(fae_.getRootCount() == fae_.connectionGraph.data.size());
	assert(root < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));

	const ConnectionGraph::CutpointSignature& signature =
		fae_.connectionGraph.data[root].signature;

	if (forbidden.count(root))
	{	// in the case the root is not to be folded
		return false;
	}

	bool found = false;
dis1_start:
	// save state offset
	fae_.pushStateOffset();

	fae_.updateConnectionGraph();

	for (const ConnectionGraph::CutpointInfo& cutpoint : signature)
	{
		if (cutpoint.root != root)
		{
			continue;
		}

		FA_DEBUG_AT(3, "type 1 cutpoint detected at root " << root);

		auto result = this->makeBox1Component(
				/* index of the TA to be folded */ root,
				/* the state where to fold */ fae_.getRoot(root)->getFinalState(),
				/* index of the other TA to be folded */ root,
				/* set of cutpoints with forbidden folding */ forbidden,
				/* box level to be discovered */ level,
				/* if true do not create the box if not present */ analysis,
				stateToBoxInfoMap
		);

		if (result.second)
		{
			found = true;
		}

		if (nullptr != result.first)
		{	// in the case folding was successful
			if (discoveredBoxes != nullptr)
			{
				discoveredBoxes->push_back(result.first);
			}

			goto dis1_start;
		}

		fae_.popStateOffset();
	}

	return found;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& s)
{
	for (const T& v : s)
		os << v << ",";

	return os;
}


bool Folding::discover2(
	size_t                       root,
	const std::set<size_t>&      forbidden,
	size_t			     level,
	bool                         analysis,
	StateToBoxInfoMap& 	    stateToBoxInfoMap,
	std::vector<const Box *>     *discoveredBoxes)
{
	// Preconditions
	assert(fae_.getRootCount() == fae_.connectionGraph.data.size());
	assert(root < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));

	const ConnectionGraph::CutpointSignature& sign =
		fae_.connectionGraph.data[root].signature;

	if (forbidden.count(root))
	{	// in the case the root is not to be folded
		return false;
	}

	bool found = false;
dis2_start:
	// save state offset
	fae_.pushStateOffset();

	fae_.updateConnectionGraph();

	// the signatures of all states in the automaton
	const ConnectionGraph::StateToCutpointSignatureMap& signatures =
		this->getSignatures(root);

	for (const ConnectionGraph::CutpointInfo& cutpoint : sign)
	{
		if (cutpoint.refCount < 2)
		{	// in the case the cutpoint is referenced less than twice
			continue;
		}

		std::unordered_set<size_t> candidateStates;

		for (const std::pair<size_t, ConnectionGraph::CutpointSignature>&
			stateSignaturePair : signatures)
		{	// for all states 's'
			for (const ConnectionGraph::CutpointInfo& tmp : stateSignaturePair.second)
			{	// for all cutpoints in the signature of state 's'
				if ((tmp.refCount < 2) /* the cutpoint is not referenced enough times from 's' */
					|| tmp.refInherited /* 's' is not a fork */
					|| (tmp.root != cutpoint.root)) /* the cutpoint is a different one */
				{	// go to the next cutpoint in the signature
					continue;
				}

				candidateStates.insert(stateSignaturePair.first);
			}
		}

		FA_DEBUG_AT(3, "candidate states: " << candidateStates);

		if (candidateStates.size() > 1)
		{
			auto ta = fae_.getRoot(root);

			std::unordered_set<size_t> visited;
			std::vector<size_t> stack(candidateStates.begin(), candidateStates.end());

			while (!stack.empty())
			{
				size_t state = stack.back();

				stack.pop_back();

				if (visited.count(state) > 0)
					continue;

				visited.insert(state);

				FA_DEBUG_AT(3, "processing: " << state);

				for (auto transitionIter = ta->begin(state); transitionIter != ta->end(state); ++transitionIter)
				{
					const Transition& t = *transitionIter;

					for (size_t child : t.GetChildren())
					{
						if (_MSB_TEST(child))
						{/* leaf state */
							continue;
						}

						FA_DEBUG_AT(3, "erasing: " << child);

						candidateStates.erase(child);
					}
				}
			}
		}

		for (size_t state : candidateStates)
		{
			FA_DEBUG_AT(3, "type 2 cutpoint detected inside component " << root
					<< " at state q" << state);

			auto result = this->makeBox1Component(
					/* index of the TA to be folded */ root,
					/* the state where to fold */ state,
					/* index of the other TA to be folded */ cutpoint.root,
					/* set of cutpoints with forbidden folding */ forbidden,
					/* box level to be discovered */ level,
					/* if true do not create the box if not present */ analysis,
					stateToBoxInfoMap
			);

			if (result.second)
			{
				found = true;
			}

			if (nullptr != result.first)
			{	// in the case folding was successful
				if (discoveredBoxes != nullptr)
				{
					discoveredBoxes->push_back(result.first);
				}

				goto dis2_start;
			}

			fae_.popStateOffset();
		}
	}

	return found;
}


bool Folding::discover3(
	size_t                      root,
	const std::set<size_t>&     forbidden,
	size_t			    level,
	bool                        analysis,
	StateToBoxInfoMap& 	    stateToBoxInfoMap,
	std::vector<const Box *>    *discoveredBoxes)
{
	// Preconditions
	assert(fae_.getRootCount() == fae_.connectionGraph.data.size());
	assert(root < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));

	const ConnectionGraph::CutpointSignature& signature =
		fae_.connectionGraph.data[root].signature;

	if (forbidden.count(root))
	{	// in the case the root is not to be folded
		return false;
	}

	bool found = false;
dis3_start:
	// save state offset
	fae_.pushStateOffset();

	fae_.updateConnectionGraph();

	for (const ConnectionGraph::CutpointInfo& cutpoint : signature)
	{	// now we search for the other cutpoint to be folded
		if (forbidden.count(cutpoint.root))
		{	// in the case the other cutpoint is not allowed to be folded
			continue;
		}

		if (root == cutpoint.root)
		{
			continue;
		}

		// retrieve the selector from 'cutpoint.root' that reaches 'root'
		size_t selectorToRoot = ConnectionGraph::getSelectorToTarget(
			fae_.connectionGraph.data[cutpoint.root].signature, root
		);

		if (selectorToRoot == static_cast<size_t>(-1))
		{	// in the case 'root' is not reachable from 'cutpoint.root'
			continue;
		}

		assert(!cutpoint.fwdSelectors.empty());

		auto result = this->makeBox2Components(
				/* index of the TA to be folded */ cutpoint.root,
				/* index of the other TA to be folded */ root,
				/* set of cutpoints with forbidden folding */ forbidden,
				/* box level to be discovered */ level,
				analysis,
				stateToBoxInfoMap,
				/* only testing that we can make the box? */ true);

		if (result.second && result.first != nullptr && !result.first->isSymetric())
		{	// in the case the box can be created in the reverse way
			FA_DEBUG_AT(3, "skipping reversed type 3 cutpoint at roots " << root << " and "
				<< cutpoint.root);
			continue;
		}

		FA_DEBUG_AT(3, "type 3 cutpoint detected at roots " << root << " and "
			<< cutpoint.root);

		result = this->makeBox2Components(
				/* index of the TA to be folded */ root,
				/* index of the other TA to be folded */ cutpoint.root,
				/* set of cutpoints with forbidden folding */ forbidden,
				/* box level to be discovered */ level,
				/* if true do not create the box if not present */ analysis,
				stateToBoxInfoMap
		);

		if (result.second)
		{
			found = true;
		}

		if (nullptr != result.first)
		{	// in the case folding was successful
			if (discoveredBoxes != nullptr)
			{
				discoveredBoxes->push_back(result.first);
			}

			goto dis3_start;
		}

		fae_.popStateOffset();
	}

	return found;
}

void Folding::componentCut(
	TreeAut&                                 res,
	TreeAut&                                 complement,
	ConnectionGraph::CutpointSignature&      complementSignature,
	const size_t                             root,
	const size_t                             state,
	const size_t                             target)
{
	// Preconditions
	assert(root < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));

	const TreeAut& src = *fae_.getRoot(root);

	res.addFinalStates(src.getFinalStates());

	complement.addFinalState(state);

	// signatures of states in the tree aut at index 'root'
	const ConnectionGraph::StateToCutpointSignatureMap& signatures = this->getSignatures(root);

	// a set of boxes that leave the state where we wish to perform the cut and
	// contain the target cutpoint in their signature
	std::unordered_set<const AbstractBox*> boxes;

	// first, we enumerate all boxes that we wish to fold
	for (auto i = src.begin(state); i != src.end(state); ++i)
	{	// traverse all transitions from 'state'
		const Transition& trans = *i;

		size_t lhsOffset = 0;

		for (const AbstractBox* box : TreeAut::GetSymbol(trans)->getNode())
		{	// go over all boxes in the transition

			// look for target cutpoint
			for (size_t j = 0; j < box->getArity(); ++j)
			{	// try all states outgoing from the box
				assert(lhsOffset + j < trans.GetChildrenSize());

				if (ConnectionGraph::containsCutpoint(
					Folding::getSignature(trans.GetNthChildren(lhsOffset + j), signatures), target))
				{	// in case the signature of the box contains the searched cutpoint
					boxes.insert(box);
					break;
				}
			}

			lhsOffset += box->getArity();
		}
	}

	for (const Transition& trans : src)
	{	// now traverse all transitions in the source tree aut
		if (trans.GetParent() != state)
		{	// if the transition does not leave 'state', simply copy it
			res.addTransition(trans);
			complement.addTransition(trans);

			continue;
		}

		// otherwise, i.e. the transition leaves 'state'

		// new LHS and label for the transition in the source automaton
		std::vector<size_t> lhs;
		std::vector<const AbstractBox*> label;

		// new LHS and label for the transition in the complement automaton
		std::vector<size_t> cLhs;
		std::vector<const AbstractBox*> cLabel;

		ConnectionGraph::CutpointSignature tmp;

		size_t lhsOffset = 0;

		// split transition
		for (const AbstractBox* box : TreeAut::GetSymbol(trans)->getNode())
		{	// traverse all boxes in the transition
			if (boxes.count(box) == 0)
			{	// in case this box does not lead to the target cutpoint, just copy it
				Folding::copyBox(lhs, label, box, trans.GetChildren(), lhsOffset);
			}
			else
			{	// in case the box _leads_ to the target cutpoint
				for (size_t j = 0; j < box->getArity(); ++j)
				{	// for each output of the box
					assert(lhsOffset + j < trans.GetChildrenSize());

					ConnectionGraph::processStateSignature(
						tmp,
						static_cast<const StructuralBox*>(box),
						j,
						trans.GetNthChildren(lhsOffset + j),
						Folding::getSignature(trans.GetNthChildren(lhsOffset + j), signatures)
					);
				}

				// copy the box into the complement
				Folding::copyBox(cLhs, cLabel, box, trans.GetChildren(), lhsOffset);
			}

			lhsOffset += box->getArity();
		}

		ConnectionGraph::normalizeSignature(tmp);

		assert(tmp.size());

		// did we arrive here for the first time?
		if (complementSignature.empty())
		{
			complementSignature = tmp;
		}

		// a bit hacky but who cares
		assert(isSignaturesCompatible(complementSignature, tmp));

		for (size_t i = 0; i < tmp.size(); ++i)
		{
			complementSignature[i].refCount =
				std::max(complementSignature[i].refCount, tmp[i].refCount);

			complementSignature[i].fwdSelectors.insert(
				tmp[i].fwdSelectors.begin(), tmp[i].fwdSelectors.end()
			);
		}

		// add the new transition to the source automaton
		assert(label.size());
		FAE::reorderBoxes(label, lhs);
		res.addTransition(lhs, fae_.boxMan->lookupLabel(label), state);

		// add the new transition to the complement automaton
		assert(cLabel.size());
		FAE::reorderBoxes(cLabel, cLhs);
		complement.addTransition(cLhs, fae_.boxMan->lookupLabel(cLabel), state);
	}
}

size_t Folding::getNestingLevel(const TreeAut& ta)
{
	size_t level = 0;
	for (auto t = ta.begin(); t != ta.end(); ++t)
	{
		auto label = TreeAut::GetSymbol(*t);
		if (!label->isNode())
		{
			continue;
		}

		for (auto b = label->getNode().begin(); b != label->getNode().end(); ++b)
			level = std::max(level, (*b)->getLevel());
	}

	return level;
}


std::pair<const Box*, bool> Folding::makeBox1Component(
		size_t root,
		size_t state,
		size_t aux,
		const std::set<size_t> &forbidden,
		size_t level,
		bool analysis,
		StateToBoxInfoMap& stateToBoxInfoMap,
		bool test)
{
	// Preconditions
	assert(root < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));

	// 'index' maintains for each cutpoint of the FA either '-1' which means that
	// the box does not reference it, or the order in which it is referenced in
	// the box
	std::vector<size_t> index(fae_.getRootCount(), static_cast<size_t>(-1));

	// maps cutpoints to selectors
	std::unordered_map<size_t, size_t> selectorMap;
	ConnectionGraph::CutpointSignature outputSignature;

	size_t start = 0;

	// split the given tree automaton at the desired location, obtain the pair of
	// the residuum and the kernel
	std::pair<TreeAutShPtr, TreeAutShPtr> resKerPair =
		this->separateCutpoint(outputSignature, root, state, aux);

	size_t boxLevel = Folding::getNestingLevel(*resKerPair.second) + 1;

	if (level != 0 && boxLevel != level)
		return std::make_pair(nullptr, boxLevel > level);

	index[root] = start++;

	for (const ConnectionGraph::CutpointInfo& cutpoint : outputSignature)
	{	// for all cutpoints in the signature of the split automaton
		if (forbidden.count(cutpoint.root))
		{	// in the case the cutpoint is forbidden to be folded
			FA_DEBUG_AT(3, "forbidden");

			return std::make_pair(nullptr, false);         // stop the folding procedure
		}

		assert(cutpoint.root < index.size());

		if (cutpoint.root != root)
		{	// in the case the cutpoint is other than the root

			// check that the signature does not contain duplicit references
			assert(static_cast<size_t>(-1) == index[cutpoint.root]);

			index[cutpoint.root] = start++;  // set the order in which it is referenced
		}
	}

	if (!Folding::computeSelectorMap(selectorMap, root, state))
	{	// in the case the box cannot be created (not all transitions from 'state'
		// have the same signature)
		FA_DEBUG_AT(3, "inconsistent signature");

		return std::make_pair(nullptr, false);
	}

	// get the input mapping of components to selector offsets
	std::vector<size_t> inputMap = extractInputMap(selectorMap, root, index);

	// create a box with a single TA
	std::unique_ptr<Box> box = std::unique_ptr<Box>(
		BoxMan::createType1Box(
			/* index of the TA put in the box */ root,
			/* the TA */ this->relabelReferences(*resKerPair.second, index),
			/* signature of the TA */ outputSignature,
			/* mapping of cutpoints to selectos */ inputMap,
			/* index renaming cutpoints */ index,
			/* nesting level */ boxLevel
		)
	);

	if (test)
	{	// in the case we are only testing
		auto boxPtr = boxMan_.lookupBox(*box);

		return std::make_pair(boxPtr, boxPtr != nullptr);
	}

	// find the box in the database
	const Box* boxPtr = this->getBox(*box, false);
/*
	if (analysis)
	{
		auto& boxInfo = stateToBoxInfoMap.insert(std::make_pair(state, std::unordered_map<const Box*, bool>())).first->second;

		assert(boxInfo.find(boxPtr) == boxInfo.end());

		boxInfo.insert(std::make_pair(boxPtr, false));

		Folding::analyzeBoxes(stateToBoxInfoMap, *box->getOutput());
	}
	else
	{
		auto stateToBoxInfoIter = stateToBoxInfoMap.find(state);

		if (stateToBoxInfoIter == stateToBoxInfoMap.end())
			return std::make_pair(nullptr, false);

		auto boxInfoIter = stateToBoxInfoIter->second.find(boxPtr);

		if (boxInfoIter == stateToBoxInfoIter->second.end() || boxInfoIter->second)
			return std::make_pair(nullptr, false);
	}
*/
	assert(boxPtr);
	assert(box->getLevel() > 1);

	FA_DEBUG_AT(2, *static_cast<const AbstractBox*>(boxPtr) << " found");

	// insert the box into the tree automaton
	fae_.setRoot(root,
		this->joinBox(
			/* the TA into which the box will be inserted */ *resKerPair.first,
			/* the state under which the box will go */ state,
			/* index of the first component in the box */ root,
			/* the box */ boxPtr,
			/* signature of the box */ outputSignature)
		);

	fae_.connectionGraph.invalidate(root);

	this->invalidateSignatures(root);

	return std::make_pair(boxPtr, true);
}


std::pair<const Box*, bool> Folding::makeBox2Components(
		size_t root,
		size_t aux,
		const std::set<size_t> &forbidden,
		size_t level,
		bool analysis,
		StateToBoxInfoMap& stateToBoxInfoMap,
		bool test)
{
	// Preconditions
	assert(root < fae_.getRootCount());
	assert(aux < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));
	assert(nullptr != fae_.getRoot(aux));

	size_t finalState = fae_.getRoot(root)->getFinalState();

	std::vector<size_t> index(fae_.getRootCount(), static_cast<size_t>(-1)), index2;
	std::vector<bool> rootMask(fae_.getRootCount(), false);

	// maps cutpoints to selectors
	std::unordered_map<size_t, size_t> selectorMap;
	ConnectionGraph::CutpointSignature outputSignature, inputSignature, tmpSignature;

	size_t start = 0;

	// split the given tree automaton at the desired location, obtain the pair of
	// the residuum and the kernel
	std::pair<TreeAutShPtr, TreeAutShPtr> resKerPair =
		this->separateCutpoint(outputSignature, root, finalState, aux);

	index[root] = start++;

	for (const ConnectionGraph::CutpointInfo& cutpoint : outputSignature)
	{
		if (cutpoint.root == root)
		{	// if this procedure cannot fold the sub-structure
			FA_DEBUG_AT(3, "refusing self-reference");

			return std::make_pair(nullptr, false);
		}

		if (forbidden.count(cutpoint.root))
		{	// in the case the cutpoint is forbidden to be folded
			FA_DEBUG_AT(3, "forbidden");

			return std::make_pair(nullptr, false);         // stop the folding procedure
		}

		assert(cutpoint.root < index.size());

		if (cutpoint.root != root)
		{	// in the case the cutpoint is other than the root

			// check that the signature does not contain duplicit references
			assert(static_cast<size_t>(-1) == index[cutpoint.root]);

			index[cutpoint.root] = start++;
		}
	}

	if (!Folding::computeSelectorMap(selectorMap, root, finalState))
	{
		FA_DEBUG_AT(3, "inconsistent signature");

		return std::make_pair(nullptr, false);
	}

	std::vector<size_t> inputMap = extractInputMap(selectorMap, root, index);

	auto auxP = this->separateCutpoint(
		inputSignature, aux, fae_.getRoot(aux)->getFinalState(), root
	);

	size_t boxLevel = std::max(
		Folding::getNestingLevel(*resKerPair.second), Folding::getNestingLevel(*auxP.second)
	) + 1;

	if (level != 0 && boxLevel != level)
		return std::make_pair(nullptr, boxLevel > level);

/*
	if (Folding::isSingular(*auxP.first))
		return false;
*/

	index2 = index;

	for (const ConnectionGraph::CutpointInfo& cutpoint : inputSignature)
	{
/*		if (cutpoint.refCount > 1)
		{
			FA_DEBUG_AT(3, "multiple references");

			return std::make_pair(nullptr, false);
		}*/

		if (forbidden.count(cutpoint.root))
		{
			FA_DEBUG_AT(3, "forbidden");

			return std::make_pair(nullptr, false);
		}

		assert(cutpoint.root < index.size());

		if (static_cast<size_t>(-1) == index[cutpoint.root])
		{
			assert(static_cast<size_t>(-1) == index2[cutpoint.root]);

			index2[cutpoint.root] = start++;

			tmpSignature.push_back(cutpoint);

			inputMap.push_back(static_cast<size_t>(-1));
		}
	}

	selectorMap.clear();

	if (!Folding::computeSelectorMap(selectorMap, aux,
		fae_.getRoot(aux)->getFinalState()))
	{
		assert(false);           // fail gracefully
	}

	size_t selector = extractSelector(selectorMap, root);

	std::shared_ptr<Box> box = std::shared_ptr<Box>(
		BoxMan::createType2Box(
			root,
			this->relabelReferences(*resKerPair.second, index),
			outputSignature,
			inputMap,
			aux,
			this->relabelReferences(*auxP.second, index2),
			inputSignature,
			selector,
			index,
			boxLevel
		)
	);

	if (test)
	{
		auto boxPtr = boxMan_.lookupBox(*box);

		return std::make_pair(boxPtr, boxPtr != nullptr);
	}

	auto boxPtr = this->getBox(*box, false);
/*
	if (analysis)
	{
		auto& boxInfo = stateToBoxInfoMap.insert(std::make_pair(finalState, std::unordered_map<const Box*, bool>())).first->second;

		assert(boxInfo.find(boxPtr) == boxInfo.end());

		boxInfo.insert(std::make_pair(boxPtr, false));

		Folding::analyzeBoxes(stateToBoxInfoMap, *box->getOutput());
		Folding::analyzeBoxes(stateToBoxInfoMap, *box->getInput());
	}
	else
	{
		auto stateToBoxInfoIter = stateToBoxInfoMap.find(finalState);

		if (stateToBoxInfoIter == stateToBoxInfoMap.end())
			return std::make_pair(nullptr, false);

		auto boxInfoIter = stateToBoxInfoIter->second.find(boxPtr);

		if (boxInfoIter == stateToBoxInfoIter->second.end() || boxInfoIter->second)
			return std::make_pair(nullptr, false);
	}
*/
	assert(boxPtr);
	assert(box->getLevel() > 1);

	FA_DEBUG_AT(2, *static_cast<const AbstractBox*>(boxPtr) << " found");

	for (const ConnectionGraph::CutpointInfo& cutpoint : tmpSignature)
	{
		outputSignature.push_back(cutpoint);
	}

	fae_.setRoot(root,
		this->joinBox(
			*resKerPair.first,
			finalState,
			root,
			boxPtr,
			outputSignature)
		);

	fae_.connectionGraph.invalidate(root);

	this->invalidateSignatures(root);

	fae_.setRoot(aux, auxP.first);
	fae_.connectionGraph.invalidate(aux);

	this->invalidateSignatures(aux);

	return std::make_pair(boxPtr, true);
}


Folding::TreeAutShPtr Folding::joinBox(
	const TreeAut&                               src,
	size_t                                       state,
	size_t                                       root,
	const Box*                                   box,
	const ConnectionGraph::CutpointSignature&    signature)
{
	// allocate a new TA and insert final states
	TreeAutShPtr ta = TreeAutShPtr(fae_.allocTA());
	ta->addFinalStates(src.getFinalStates());

	for (const Transition& trans : src)
	{	// for every transition in the source
		if (trans.GetParent() != state)
		{	// in the case the state is not important, simply copy the transition
			ta->addTransition(trans);

			continue;
		}

		std::vector<const AbstractBox*> label(TreeAut::GetSymbol(trans)->getNode());
		std::vector<size_t> lhs(trans.GetChildren());

		// append the box
		label.push_back(box);

		for (const ConnectionGraph::CutpointInfo& cutpoint : signature)
		{	// for all cutpoints in the signature (with the exception of those that
			// are inside the box), append references to other (redundant) cutpoints
			// the 'lhs'
			if ((cutpoint.root == root) && (src.getFinalStates().count(state)))
			{	// in the case the cutpoint is the box in the joined automaton or
				// 'state' is final
				continue;
			}

			lhs.push_back(
				fae_.addData(*ta, Data::createRef(cutpoint.root))
			);
		}

		// reorder the boxes in 'label' and states in 'lhs' so that it is ordered 
		FA::reorderBoxes(label, lhs);

		ta->addTransition(lhs, fae_.boxMan->lookupLabel(label), state);
	}

	return ta;
}


bool Folding::checkSelectorMap(
	const std::unordered_map<size_t, size_t>&     selectorMap,
	size_t                                        root,
	size_t                                        state)
{
	// Preconditions
	assert(root < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));

	const ConnectionGraph::StateToCutpointSignatureMap& signatures =
		this->getSignatures(root);

	const TreeAut& ta = *fae_.getRoot(root);

	for (auto i = ta.begin(state); i != ta.end(state); ++i)
	{
		std::unordered_map<size_t, size_t> m;

		computeSelectorMapTrans(m, *i, signatures);

		if (m != selectorMap)
		{
			return false;
		}
	}

	return true;
}


bool Folding::computeSelectorMap(
	std::unordered_map<size_t, size_t>&      selectorMap,
	const size_t                             root,
	const size_t                             state)
{
	// Preconditions
	assert(root < fae_.getRootCount());
	assert(nullptr != fae_.getRoot(root));
	assert(selectorMap.empty());

	const TreeAut& ta = *fae_.getRoot(root);

	// check there are some transitions from 'state'
	assert(ta.begin(state) != ta.end(state));

	const ConnectionGraph::StateToCutpointSignatureMap& signatures =
		this->getSignatures(root);

	computeSelectorMapTrans(selectorMap, *ta.begin(state), signatures);

	// check whether the computed selector map is the same for all transitions
	// leaving 'state'
	return this->checkSelectorMap(selectorMap, root, state);
}

void Folding::learn1(FAE& fae, BoxMan& boxMan, std::set<size_t> forbidden)
{
	throw "not implemented";
}

void Folding::learn2(FAE& fae, BoxMan& boxMan, std::set<size_t> forbidden)
{
	throw "not implemented";
}

std::pair<std::unordered_map<size_t, std::vector<const Box *>>, bool> Folding::fold(
        FAE&                         fae,
        BoxMan&                      boxMan,
#if FA_BOX_APPROXIMATION
	BoxAntichain&		     boxAntichain,
#endif
        const std::set<size_t>&      forbidden,
        size_t			     level,
        bool			     discover3Only,
        bool			     analysis,
        StateToBoxInfoMap&	     stateToBoxInfoMap
        )
{
	std::unordered_map<size_t, std::vector<const Box *>> foldedRoots;

#if FA_BOX_APPROXIMATION
//	Folding folding(fae, boxMan, boxMan.boxAntichain());
	Folding folding(fae, boxMan, boxAntichain);
#else
	Folding folding(fae, boxMan);
#endif

	bool found = false;

	for (size_t i = 0; i < fae.getRootCount(); ++i)
	{
		if (forbidden.end() != forbidden.find(i))
		{    // in the case the cutpoint is not allowed for folding
			continue;
		}

		assert(nullptr != fae.getRoot(i));

		// Try to fold the 3 types of cutpoints starting from cutpoint 'i', but
		// _ONLY_ using boxes which are _ALREADY_ in 'boxMan'. No learning of new
		// boxes is allowed

		foldedRoots[i] = std::vector<const Box *>();
		found |= folding.discover3(i, forbidden, level, analysis, stateToBoxInfoMap, &foldedRoots.at(i));
		if (!discover3Only)
		{
			found |= folding.discover1(i, forbidden, level, analysis, stateToBoxInfoMap, &foldedRoots.at(i));
			found |= folding.discover2(i, forbidden, level, analysis, stateToBoxInfoMap, &foldedRoots.at(i));
		}

		if (foldedRoots.at(i).size() == 0)
		{
			foldedRoots.erase(i);
		}
	}

	if (foldedRoots.size())
	{
		FA_DEBUG_AT(3, "after folding: " << std::endl << fae);
	}

	return std::make_pair(foldedRoots, found);
}

bool Folding::updateBoxSignature(
		StateToBoxSignatureMap& stateMap,
		size_t             		state,
		const BoxSignature&		boxSignature)
{
	auto itBoolPair = stateMap.insert(std::make_pair(state, boxSignature));

	if (itBoolPair.second)
		return true;

	// in case the state is already mapped to something
	BoxSignature& prevSignature = itBoolPair.first->second;

	bool changed = false;

	for (auto& boxOriginPair : prevSignature)
	{
		auto iter = boxSignature.find(boxOriginPair.first);
		if (iter == boxSignature.end() && (boxOriginPair.second != static_cast<size_t>(-1)))
		{
			boxOriginPair.second = static_cast<size_t>(-1);
			changed = true;
		}
	}

	for (auto& boxOriginPair : boxSignature)
	{
		auto itBoolPair2 = prevSignature.insert(boxOriginPair);

		if (itBoolPair2.second || (boxOriginPair.second != itBoolPair2.first->second && itBoolPair2.first->second != static_cast<size_t>(-1)))
		{
			itBoolPair2.first->second = static_cast<size_t>(-1);
			changed = true;
		}
	}

	return changed;
}


void Folding::joinBoxSignature(BoxSignature& dst, const std::pair<const Box*, size_t>& boxOriginPair)
{
	auto itBoolPair = dst.insert(boxOriginPair);

	if (itBoolPair.second)
		return;

	// -1 indicates that box appears multiple times
	itBoolPair.first->second = static_cast<size_t>(-1);
}


void Folding::joinBoxSignature(BoxSignature& dst, const BoxSignature& src)
{
	for (auto& boxOriginPair : src)
		joinBoxSignature(dst, boxOriginPair);
}


bool Folding::processLhs(
		BoxSignature&					result,
		const std::vector<size_t>&  	lhs,
		const StateToBoxSignatureMap&	stateMap)
{
	for (auto state : lhs)
	{
		auto it = stateMap.find(state);

		if (it == stateMap.end())
			return false;

		joinBoxSignature(result, it->second);
	}

	return true;
}


void Folding::computeBoxSignatureMap(
		StateToBoxSignatureMap&	stateMap,
		const TreeAut&			ta)
{
	stateMap.clear();

	// the workset of transitions
	std::vector<TreeAut::Transition> transitions;

	// compute the initial signatures for leaves, other signatures are cleared
	for (auto trans : ta)
	{	// traverse transitions of the TA
		if (TreeAut::GetSymbol(trans)->isData())
		{	// for data transitions
			stateMap.insert(std::make_pair(trans.GetParent(), BoxSignature()));
		} else
		{	// for non-data transitions
			transitions.push_back(trans);
		}
	}

	BoxSignature boxSignature;

	// Now we propagate the computed signatures upward in the tree structure until the signatures
	// stabilize. 'transitions' contains transitions that are to be processed.
	bool changed = true;
	while (changed)
	{	// while there are still some transitions to be processed
		changed = false;
		for (auto t : transitions)
		{
			boxSignature.clear();

			if (!processLhs(boxSignature, t.GetChildren(), stateMap))
			{	// in case this transition cannot be processed because of some downward
				// states with missing box information
				continue;
			}

			for (const AbstractBox* box : TreeAut::GetSymbol(t)->getNode())
			{	// for all boxes in the label
				assert(nullptr != box);

				if (box->isBox())
				{
					joinBoxSignature(boxSignature, std::make_pair(static_cast<const Box*>(box), t.GetParent()));
				}
			}

			changed |= updateBoxSignature(stateMap, t.GetParent(), boxSignature);
		}
	}
}

void Folding::analyzeBoxes(StateToBoxInfoMap& stateToBoxInfoMap, const TreeAut& ta)
{
	StateToBoxSignatureMap stateMap;

	computeBoxSignatureMap(stateMap, ta);

	for (auto finalState : ta.getFinalStates())
	{
		auto boxSignatureIter = stateMap.find(finalState);

		assert(boxSignatureIter != stateMap.end());

		for (auto& boxStatePair : boxSignatureIter->second)
		{
			if (boxStatePair.second == static_cast<size_t>(-1))
				continue;

			auto stateBoxInfoIter = stateToBoxInfoMap.find(boxStatePair.second);

			if (stateBoxInfoIter == stateToBoxInfoMap.end())
				continue;

			auto boxInfoIter = stateBoxInfoIter->second.find(boxStatePair.first);

			if (boxInfoIter == stateBoxInfoIter->second.end())
				continue;

			boxInfoIter->second = true;
		}
	}
}
