/*
 * Copyright (C) 2012 Ondrej Lengal, Martin Hruska
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

/**
 * @file symstate_isect.cc
 *
 * This file contains implementation of functions of SymState that perform
 * intersection of a pair of symbolic states.
 */

#include <assert.h>

// Forester headers
#include "streams.hh"
#include "symstate.hh"
#include "virtualmachine.hh"

namespace
{
	std::string stateToString(const size_t state)
	{
		std::ostringstream os;

		if (_MSB_TEST(state))
			os << 'r' << _MSB_GET(state);
		else
			os << 'q' << state;

		return os.str();
	}
}

struct RootState
{
	size_t root;
	size_t state;

	RootState() :
		root(0),
		state(0)
	{ }

	RootState(size_t pRoot, size_t pState) :
		root(pRoot),
		state(pState)
	{ }

	bool operator<(const RootState& rhs) const
	{
		if (root < rhs.root)
		{
			return true;
		}
		else if (root == rhs.root)
		{
			return state < rhs.state;
		}
		else
		{
			return false;
		}
	}

	friend std::ostream& operator<<(std::ostream& os, const RootState& rs)
	{
		os << "[" << rs.root << ", " << FA::writeState(rs.state) << "]";
		return os;
	}
};


typedef std::pair<RootState, RootState> ProdState;


// anonymous namespace
namespace
{
/**
 * @brief  Engine for reference substitution
 */
class SubsRefEngine
{
private:  // data members

	/// Keeps product states, e.g. (p,q)
	std::set<ProdState> processed_;

	/// The work stack
	std::vector<ProdState> workstack_;

public:   // methods

	SubsRefEngine() :
		processed_(),
		workstack_()
	{ }

	/**
	 * @brief  Creates a product state and adds it to processing
	 *
	 * Creates a product state in case it does not exist.
	 *
	 * @param[in]  lhsRoot   Root number of the LHS automaton
	 * @param[in]  lhsState  State in the LHS automaton
	 * @param[in]  rhsRoot   Root number of the RHS automaton
	 * @param[in]  rhsState  State in the RHS automaton
	 */
	void makeProductState(
		const size_t&          lhsRoot,
		const size_t&          lhsState,
		const size_t&          rhsRoot,
		const size_t&          rhsState,
		const bool&            jumped=false)
	{
		// the product state
		RootState lhsRootState(lhsRoot, lhsState);
		RootState rhsRootState(rhsRoot, rhsState);
		ProdState prodState(lhsRootState, rhsRootState);

		auto itBoolPairProcessed = processed_.insert(prodState);
		bool isNewState = itBoolPairProcessed.second;

		if (isNewState)
		{	// in case the state has not been processed before
			FA_DEBUG_AT(1,"Creating new product state (" <<  lhsRootState << ", "
				<< rhsRootState << ")");
			workstack_.push_back(prodState);
		}
	}


	/**
	 * @brief  Determines whether the work stack is empty
	 *
	 * This method returns @p true in case the work stack is empty, @p false
	 * otherwise.
	 *
	 * @returns  @p true in case the work stack is empty, @p false otherwise
	 */
	bool wsEmpty() const
	{
		return workstack_.empty();
	}

	/**
	 * @brief  Retrieves the next state to be processed
	 *
	 * This method retrieves the next state that has not been processed so far.
	 * Further, it removes it from the work stack.
	 *
	 * @returns  The product state
	 */
	ProdState getNextState()
	{
		// Preconditions
		assert(!wsEmpty());

		ProdState res = *workstack_.crbegin();
		workstack_.pop_back();
		return res;
	}

	const std::set<ProdState>& getProcessed() const
	{
		return processed_;
	}
};

/**
 * @brief  Engine for intersection
 */
class IsectEngine
{
private:  // data members

	/// FAE to be the output
	FAE& fae_;

	/// Maps product states to states in the new automaton, e.g. (p,q) -> r
	std::map<ProdState, RootState> processed_;

	/// The work stack
	std::vector<std::pair<ProdState, RootState>> workstack_;

	/// Maps pairs of roots to a root in the new automaton
	std::map<std::pair<size_t, size_t>, size_t> rootMap_;

	/// Maps roots of RHS automaton to roots of the new automaton
	std::map<size_t, size_t> rhsRootMap_;

	/// Counter of roots in the new FAE
	size_t rootCnt_;

public:   // methods

	explicit IsectEngine(FAE& fae) :
		fae_(fae),
		processed_(),
		workstack_(),
		rootMap_(),
		rhsRootMap_(),
		rootCnt_(0)
	{ }

	/**
	 * @brief  Adds a state to the automaton if it does not exist
	 *
	 * Inserts a product state in the automaton in case it does not exist. In
	 * either case returns the number of the root and the number of the state in
	 * the target FAE (in the form of a pair). In case it is the first state of
	 * a new automaton (previously unknown combination of root numbers), it is set
	 * as the final state.
	 *
	 * @param[in]  lhsRoot   Root number of the LHS automaton
	 * @param[in]  lhsState  State in the LHS automaton
	 * @param[in]  rhsRoot   Root number of the RHS automaton
	 * @param[in]  rhsState  State in the RHS automaton
	 *
	 * @returns  Root in the new FAE and the number of the product state
	 */
	RootState makeProductState(
		const size_t&          lhsRoot,
		const size_t&          lhsState,
		const size_t&          rhsRoot,
		const size_t&          rhsState,
		const bool&            jumped=false)
	{
		// the product state
		RootState lhsRootState(lhsRoot, lhsState);
		RootState rhsRootState(rhsRoot, rhsState);
		ProdState prodState(lhsRootState, rhsRootState);

		// find the root # in the aut (or create one)
		auto itBoolPairRootMap = rootMap_.insert(std::make_pair(
			std::make_pair(lhsRoot, rhsRoot), rootCnt_));
		bool isNewRoot = itBoolPairRootMap.second;
		if (isNewRoot)
		{
			FA_DEBUG_AT(1,"Creating new root: " << rootCnt_ << " as the product of roots ("
				<< lhsRoot << ", " << rhsRoot << ")");
			++rootCnt_;

			TreeAut* ta = fae_.allocTA();
			fae_.appendRoot(ta);
			fae_.connectionGraph.newRoot();
		}

		// the actual number of the root
		const size_t& root = itBoolPairRootMap.first->second;

		RootState newState(root, fae_.nextState());
		auto itBoolPairProcessed = processed_.insert(std::make_pair(prodState, newState));
		bool isNewState = itBoolPairProcessed.second;

		FA_DEBUG_AT(1,"Processing product state [" << rhsState << "," << lhsState << "] -> " 
				<< newState << ", roots ["
				<< lhsRoot	<< "," << rhsRoot << "] -> " << root);

		// isNewRoot -> isNewState
		assert(!isNewRoot || isNewState);

		if (isNewState)
		{	// in case the state has not been processed before
			FA_DEBUG_AT(1,"Creating new state: " << FA::writeState(fae_.nextState())
				<< " as the product of states (" <<  lhsRootState << ", "
				<< rhsRootState << ")");
			workstack_.push_back(std::make_pair(prodState, newState));
			fae_.newState();
		}

		// the actual number of the state
		const size_t& state = itBoolPairProcessed.first->second.state;

		if (isNewRoot || jumped)
		{	// set final state
			FA_DEBUG_AT(1,"Creating new final state " << state << "("
				<< itBoolPairProcessed.first->first.first << ","
				<< itBoolPairProcessed.first->first.second << ")"
				<< " of root " << root); 
			fae_.getRoot(root)->addFinalState(state);
		}

		auto itBoolRhsRootMap = rhsRootMap_.insert(std::make_pair(rhsRoot, root));
		if (!itBoolRhsRootMap.second)
		{
			FA_DEBUG_AT(1,"Not entering already processed RHS root " << rhsRoot << " -> "
				<< root << " (already present mapping "
				<< itBoolRhsRootMap.first->first << " -> "
				<< itBoolRhsRootMap.first->second << ")");
		}

		return RootState(root, state);
	}


	/**
	 * @brief  Determines whether the work stack is empty
	 *
	 * This method returns @p true in case the work stack is empty, @p false
	 * otherwise.
	 *
	 * @returns  @p true in case the work stack is empty, @p false otherwise
	 */
	bool wsEmpty() const
	{
		return workstack_.empty();
	}

	/**
	 * @brief  Retrieves the next state to be processed
	 *
	 * This method retrieves the next state that has not been processed so far.
	 * Further, it removes it from the work stack.
	 *
	 * @returns  A pair of a product state and corresponding root and state in new
	 *           FAE
	 */
	std::pair<ProdState, RootState> getNextState()
	{
		// Preconditions
		assert(!wsEmpty());

		std::pair<ProdState, RootState> res = *workstack_.crbegin();
		workstack_.pop_back();
		return res;
	}

	std::vector<size_t> getRootOrderIndexForRHS() const
	{
		std::vector<size_t> index(fae_.getRootCount(), static_cast<size_t>(-1));

		std::set<size_t> usedVals;

		for (const std::pair<size_t, size_t>& rootPair : rhsRootMap_)
		{
			const size_t& oldRoot = rootPair.first;
			const size_t& newRoot = rootPair.second;
			assert(static_cast<size_t>(-1) == index[newRoot]);
			index[newRoot] = oldRoot;

			if (!usedVals.insert(oldRoot).second)
			{	// in case the value could not be inserted
				assert(false);               // fail gracefully
			}
		}

		size_t cnt = 0;

		for (size_t i = 0; i < index.size(); ++i)
		{	// fill not bound root numbers
			if (static_cast<size_t>(-1) == index[i])
			{	// in case the root number has not been bound
				while (usedVals.find(cnt) != usedVals.cend())
				{	// try to find a non-used root number
					++cnt;
				}

				FA_DEBUG_AT(1,"Adding " << cnt << " at position " << i << " of index");
				index[i] = cnt++;
			}
		}

		return index;
	}
};

bool shouldCreateProductState(
		const FAE& thisFAE,
		const FAE& srcFAE,
		const Data& thisVar,
		const Data& srcVar)
{
	auto isundef = [](const Data& data) -> bool {return data.isUndef();};
	auto isnative = [](const Data& data) -> bool {return data.isNativePtr();};

	if (VirtualMachine::isNodeType(thisFAE, thisVar, isundef) &&
		VirtualMachine::isNodeType(srcFAE, srcVar, isundef))
	{ // variable in register is undefined pointer so it has not its own TA
		return false;
	}
	else if (VirtualMachine::isNodeType(thisFAE, thisVar, isnative) &&
		VirtualMachine::isNodeType(srcFAE, srcVar, isnative))
	{ // native pointer is return address of a function
		return false;
	}

	return true;
}

} // namespace


void SymState::SubstituteRefs(
	const SymState&      src,
	const Data&          oldValue,
	const Data&          newValue)
{
	// Preconditions
	assert(oldValue.isRef() && newValue.isRef());

	const std::shared_ptr<const FAE> thisFAE = this->GetFAE();
	const std::shared_ptr<const FAE> srcFAE = src.GetFAE();
	assert((nullptr != thisFAE) && (nullptr != srcFAE));

	FA_DEBUG_AT(1,"before substitution: " << *thisFAE);

	FAE* fae = new FAE(*thisFAE);
	fae->clear();
	fae->setStateOffset(thisFAE->nextState());
	this->SetFAE(std::shared_ptr<FAE>(fae));

	// engine that handles creation of new states etc.
	SubsRefEngine engine;

	for (size_t i = 0; i < thisFAE->getRootCount(); ++i)
	{	// allocate existing TA
		if (nullptr != thisFAE->getRoot(i))
		{
			TreeAut* ta = fae->allocTA();
			ta->addFinalStates(thisFAE->getRoot(i)->getFinalStates());
			fae->appendRoot(ta);
		}
		else
		{
			fae->appendRoot(nullptr);
		}

		fae->connectionGraph.newRoot();
	}

	assert(thisFAE->GetVarCount() == srcFAE->GetVarCount());
	for (size_t i = 0; i < thisFAE->GetVarCount(); ++i)
	{	// copy global variables
		Data thisVar = thisFAE->GetVar(i);
		const Data& srcVar = srcFAE->GetVar(i);
		const bool wasThisUndef = thisVar.isUndef();

		if (srcVar.isRef() && wasThisUndef)
		{	// in case we need to substitute at global variable
			// what could happen when a node was referenced by the global
			// variable or there was a value returned by a function.
			assert(oldValue == srcVar || i > 1);
			thisVar = newValue;
		}

		fae->PushVar(thisVar);

		if (!srcVar.isRef() || !thisVar.isRef())
		{	// in case some of them is not a reference
			assert(srcVar == thisVar);
			continue;
		}

		if (wasThisUndef || i > 1) // TODO is this really correct
		{ // the value in a global register of thisFA contained a value returned by function
			continue;
		}

		assert(thisVar.isRef() && srcVar.isRef());
		assert(thisVar.d_ref.displ == srcVar.d_ref.displ);
		assert(0 == thisVar.d_ref.displ);

		const TreeAut* thisRoot = thisFAE->getRoot(thisVar.d_ref.root).get();
		const TreeAut* srcRoot  = srcFAE->getRoot(srcVar.d_ref.root).get();
		assert((nullptr != thisRoot) && (nullptr != srcRoot));

		engine.makeProductState(
			thisVar.d_ref.root, thisRoot->getFinalState(),
			srcVar.d_ref.root, srcRoot->getFinalState());
	}

	assert(this->GetRegCount() == src.GetRegCount());
	for (size_t i = 0; i < this->GetRegCount(); ++i)
	{	// copy registers
		Data thisVar = this->GetReg(i);
		const Data& srcVar = src.GetReg(i);

		if (srcVar.isRef() && thisVar.isUndef())
		{	// in case we need to substitute at global variable
			assert(oldValue == srcVar);
			this->SetReg(i, thisVar);
			thisVar = newValue;
		}

		if (!srcVar.isRef() || !thisVar.isRef())
		{	// in case some of them is not a reference
			assert(srcVar == thisVar);
			continue;
		}

		assert(thisVar.isRef() && srcVar.isRef());
		assert(thisVar.d_ref.displ == srcVar.d_ref.displ);

		size_t thisRef = thisVar.d_ref.root;
		size_t srcRef = srcVar.d_ref.root;

		if (thisVar.d_ref.displ != 0)
		{ // not the TA referenced by thisVar is used to make a new product state
			auto isref = [](const Data& data) -> bool {return data.isRef();};

			if (!shouldCreateProductState(*thisFAE, *srcFAE, thisVar, srcVar))
			{
				continue;
			}
			else if (VirtualMachine::isNodeType(*thisFAE, thisVar, isref) &&
				VirtualMachine::isNodeType(*srcFAE, srcVar, isref))
			{ // TODO is this correct?
				// if both register points a reference than the referenced root is loaded
				FA_WARN("Suspicious replacement of root refs");
				Data tmpData;
				VirtualMachine(*thisFAE).nodeLookup(thisVar.d_ref.root, thisVar.d_ref.displ, tmpData);
				thisRef = tmpData.d_ref.root;
				VirtualMachine(*srcFAE).nodeLookup(thisVar.d_ref.root, thisVar.d_ref.displ, tmpData);
				srcRef = tmpData.d_ref.root;
			}
			else
			{
				assert(false);
			}
		}

		if (!thisFAE->rootDefined(thisRef) && !srcFAE->rootDefined(srcRef))
		{ // register contains reference which has been removed before this instruction in forward run
			continue;
		}

		const TreeAut* thisRoot = thisFAE->getRoot(thisRef).get();
		const TreeAut* srcRoot  = srcFAE->getRoot(srcRef).get();

		assert((nullptr != thisRoot) && (nullptr != srcRoot));

		engine.makeProductState(
			thisVar.d_ref.root, thisRoot->getFinalState(),
			srcVar.d_ref.root, srcRoot->getFinalState());
	}

	while (!engine.wsEmpty())
	{
		const ProdState curState = engine.getNextState();

		const size_t& thisRoot = curState.first.root;
		const size_t& srcRoot = curState.second.root;

		const size_t& thisState = curState.first.state;
		const size_t& srcState = curState.second.state;

		FA_DEBUG_AT(1,"Processing product state (" <<  curState.first << ", "
			<< curState.second << ")");

		const std::shared_ptr<TreeAut> thisTA = thisFAE->getRoot(thisRoot);
		const std::shared_ptr<TreeAut> srcTA = srcFAE->getRoot(srcRoot);
		assert((nullptr != thisTA) && (nullptr != srcTA));

		auto thisIt = thisTA->begin(thisState);
		auto thisEnd = thisTA->end(thisState);
		auto srcEnd = srcTA->end(srcState);

		for (; thisIt != thisEnd; ++thisIt)
		{
			auto srcIt = srcTA->begin(srcState);
			for (; srcIt != srcEnd; ++srcIt)
			{
				const auto& thisTrans = *thisIt;
				const auto& srcTrans = *srcIt;

				// we handle data one level up
				assert(!TreeAut::GetSymbol(thisTrans)->isData() 
                        && !TreeAut::GetSymbol(srcTrans)->isData());
				assert(!thisTrans.GetChildren().empty() && !srcTrans.GetChildren().empty());

				// TODO: so far, we are not doing unfolding!
				if (TreeAut::GetSymbol(thisTrans) == TreeAut::GetSymbol(srcTrans))
				{
					FA_DEBUG_AT(1,"Transition: " << 
                            TreeAut::GetSymbol(thisTrans) <<
                            ", thisState = " <<
                            thisState <<
                            ", srcState = " << srcState);

					assert(thisTrans.GetChildrenSize() == srcTrans.GetChildrenSize());
					const size_t& transArity = thisTrans.GetChildrenSize();

					if (TreeAut::GetSymbol(thisTrans)->isData())
					{	// data are processed one level up
						assert(TreeAut::GetSymbol(srcTrans)->isData());
						continue;
					}

					std::vector<size_t> lhs;
					size_t i;
					for (i = 0; i < transArity; ++i)
					{	// for each pair of states that map to each other
						const Data* srcData = nullptr, *thisData = nullptr;
						const bool  srcIsData =
								srcFAE->isData(srcTrans.GetNthChildren(i),  srcData);
						const bool thisIsData =
								thisFAE->isData(thisTrans.GetNthChildren(i), thisData);

						//  srcIsData <-> nullptr != srcData
						assert((!srcIsData || (nullptr != srcData))
							&& (srcIsData || (nullptr == srcData)));
						//  thisIsData <-> nullptr != thisData
						assert((!thisIsData || (nullptr != thisData))
							&& (thisIsData || (nullptr == thisData)));

						if (!srcIsData && !thisIsData)
						{	// ************* process internal states *************
							// This is the easiest case, when both states in the product are
							// internal.
							engine.makeProductState(
								thisRoot, thisTrans.GetNthChildren(i),
								srcRoot, srcTrans.GetNthChildren(i));

							lhs.push_back(thisTrans.GetNthChildren(i));
						}
						else if (srcIsData && thisIsData &&
							!srcData->isRef() && !thisData->isRef())
						{ // ************* process real data states (leaves) *************
							// This is the second easiest case, when both states are real data
							// (i.e. no references). In this case, we are simply doing
							// intersection of the data.

							// Ignore if it is not compatible... we will detect it at
							// intersection
							//assert(*thisData == *srcData);

							lhs.push_back(fae->addData(
								*fae->getRoot(thisRoot).get(), *thisData));
						}
						else if (srcIsData && thisIsData &&
							srcData->isRef() && thisData->isRef())
						{ // ************* process reference states (leaves) *************
							// This is another quite easy case, when both states are
							// references to another automata. In this case, we are jumping
							// from both automata into another product automaton.
							assert(thisData->d_ref.displ == srcData->d_ref.displ);
							assert(0 == thisData->d_ref.displ);

							FA_DEBUG_AT(1,"Two references");

							const size_t& thisNewRoot = thisData->d_ref.root;
							const size_t& srcNewRoot  = srcData->d_ref.root;

							const TreeAut* thisNewTA = thisFAE->getRoot(thisNewRoot).get();
							const TreeAut* srcNewTA  = srcFAE->getRoot(srcNewRoot).get();
							assert((nullptr != thisNewTA) && (nullptr != srcNewTA));

							engine.makeProductState(
								thisNewRoot, thisNewTA->getFinalState(),
								srcNewRoot, srcNewTA->getFinalState());

							lhs.push_back(fae->addData(*fae->getRoot(thisRoot).get(),
								Data::createRef(thisNewRoot)));
						}
						else if (srcIsData && (*srcData == oldValue))
						{ // ************* perform substitution of reference *************
							assert(thisIsData && thisData->isUndef());

							FA_DEBUG_AT(1,"Substituting " << *srcData << " for " << newValue);

							lhs.push_back(fae->addData(
								*fae->getRoot(thisRoot).get(), newValue));
						}
						else if ((srcIsData && !thisIsData && srcData->isNull())
							|| (!srcIsData && thisIsData && thisData->isNull())
							|| (srcIsData && thisIsData && srcData->isNull() && thisData->isRef())
							|| (srcIsData && thisIsData && srcData->isRef() && thisData->isNull()))
						{ // ************* process NULL pointers *************
							// This is the case when there is a NULL pointer and either an
							// internal state or a reference
							FA_DEBUG_AT(1,"NULL ptr!");

							break;   // cut this branch of the product
						}
						else if ((srcIsData && !thisIsData && srcData->isRef())
							|| (!srcIsData && thisIsData && thisData->isRef()))
						{ // ************* process jumps *************
							// This is the case when one FA jumps to another TA and the other does not
							FA_DEBUG_AT(1,"jump!");

							if (srcIsData)
							{
								assert(srcIsData && !thisIsData && srcData->isRef());

								const size_t& srcNewRoot = srcData->d_ref.root;
								const TreeAut* srcNewTA  = srcFAE->getRoot(srcNewRoot).get();
								assert(nullptr != srcNewTA);

								engine.makeProductState(
									thisRoot, thisTrans.GetNthChildren(i),
									srcNewRoot, srcNewTA->getFinalState(),
									true);

								lhs.push_back(thisTrans.GetNthChildren(i));
							}
							else
							{
								assert(!srcIsData && thisIsData && thisData->isRef());

								const size_t& thisNewRoot = thisData->d_ref.root;
								const TreeAut* thisNewTA = thisFAE->getRoot(thisNewRoot).get();
								assert(nullptr != thisNewTA);

								engine.makeProductState(
									thisNewRoot, thisNewTA->getFinalState(),
									srcRoot, srcTrans.GetNthChildren(i),
									true);

								lhs.push_back(fae->addData(*fae->getRoot(thisRoot).get(),
									Data::createRef(thisNewRoot)));
							}
						}
						else
						{	// we should not get here
							assert(false);       // fail gracefully
						}
					}

					if (transArity == i)
					{	// in case we have not interrupted the search, add the transition
						std::ostringstream osLhs;
						for (auto it = lhs.cbegin(); it != lhs.cend(); ++it)
						{
							if (lhs.cbegin() != it)
								osLhs << ",";

							osLhs << FA::writeState(*it);
						}

						FA_DEBUG_AT(1,"TA " << thisRoot << ": adding transition "
							<< FA::writeState(thisState) << " -> "
							<< TreeAut::GetSymbol(thisTrans) << "(" << osLhs.str() << ")");

						fae->getRoot(thisRoot)->addTransition(lhs,
								TreeAut::GetSymbol(thisTrans),
								thisState);
					}
				}
				else
				{
					FA_DEBUG_AT(1,"not-matching: " << *thisIt << ", " << *srcIt);
				}
			}
		}
	}

	// TODO: 'odprasit' this function
	// TODO: I also need to perform some on-the fly splitting
	FA_DEBUG_AT(1,"after substitution: " << *fae);

//	fae->updateConnectionGraph();
}


void SymState::Intersect(
	const SymState&          fwd)
{
	const std::shared_ptr<const FAE> thisFAE = this->GetFAE();
	const std::shared_ptr<const FAE> fwdFAE = fwd.GetFAE();
	assert((nullptr != thisFAE) && (nullptr != fwdFAE));

	FAE* fae = new FAE(*thisFAE);
	fae->clear();
	fae->setStateOffset(thisFAE->nextState());
	this->SetFAE(std::shared_ptr<FAE>(fae));

	// engine that handles creation of new states etc.
	IsectEngine engine(*fae);

	if (this->GetRegCount() != fwd.GetRegCount())
	{	// if the number of local registers does not match
		FA_DEBUG_AT(1, "Number of local registers does not match -> creating empty intersection");
		return;      // empty FA
	}

	if (thisFAE->GetVarCount() != fwdFAE->GetVarCount())
	{	// if the number of input ports of the FAE does not match
		FA_DEBUG_AT(1, "Number of input ports does not match -> creating empty intersection");
		return;      // empty FA
	}

	assert(thisFAE->GetVarCount() == fwdFAE->GetVarCount());
	for (size_t i = 0; i < thisFAE->GetVarCount(); ++i)
	{	// check global variables
		const Data& thisVar = thisFAE->GetVar(i);
		const Data& fwdVar = fwdFAE->GetVar(i);

		if (!thisVar.isRef() && !fwdVar.isRef())
		{	// in case of non-references
			if (thisVar != fwdVar)
			{
				return;   // empty FA
			}
		}
		else
		{
			assert(thisVar.isRef() && fwdVar.isRef());
		}
	}

	assert(this->GetRegCount() == fwd.GetRegCount());
	for (size_t i = 0; i < this->GetRegCount(); ++i)
	{	// check local registers

		// NOTE: this needs to be done in order to also collect temporary TA
		// references for parts of heap which has not been so far connected so that
		// they would be reachable from global variables

		const Data& thisVar = this->GetReg(i);
		const Data& fwdVar = fwd.GetReg(i);

		if (!thisVar.isRef() && !fwdVar.isRef())
		{	// in case of non-references
			if (thisVar != fwdVar)
			{
				return;   // empty FA
			}
		}
		else
		{
			assert(thisVar.isRef() && fwdVar.isRef());
		}
	}

	assert(thisFAE->GetVarCount() == fwdFAE->GetVarCount());
	for (size_t i = 0; i < thisFAE->GetVarCount(); ++i)
	{	// add processing of all global variables
		const Data& thisVar = thisFAE->GetVar(i);
		const Data& fwdVar = fwdFAE->GetVar(i);

		if (!thisVar.isRef() || !fwdVar.isRef())
		{	// in case some of them is not a reference
			assert(fwdVar == thisVar);
			fae->PushVar(thisVar);
		}
		else if (0 != thisVar.d_ref.displ &&
				!shouldCreateProductState(*thisFAE, *fwdFAE, thisVar, fwdVar))
		{
			assert(thisVar.d_ref.displ == fwdVar.d_ref.displ);
			fae->PushVar(thisVar);
		}
		else
		{
			assert(thisVar.isRef() && fwdVar.isRef());
			assert(thisVar.d_ref.displ == fwdVar.d_ref.displ);
			assert(0 == thisVar.d_ref.displ);

			const size_t thisRootNum = thisVar.d_ref.root;
			const size_t fwdRootNum = fwdVar.d_ref.root;

			const TreeAut* thisRoot = thisFAE->getRoot(thisRootNum).get();
			const TreeAut* fwdRoot  = fwdFAE->getRoot(fwdRootNum).get();
			assert((nullptr != thisRoot) && (nullptr != fwdRoot));

			RootState rs = engine.makeProductState(
				thisRootNum, thisRoot->getFinalState(),
				fwdRootNum, fwdRoot->getFinalState()
			);

			fae->PushVar(Data::createRef(rs.root));
		}
	}

	while (!engine.wsEmpty())
	{
		const std::pair<ProdState, RootState> curProdState = engine.getNextState();

		const ProdState& curState = curProdState.first;
		const RootState& curNewState = curProdState.second;

		const size_t& thisRoot = curState.first.root;
		const size_t& fwdRoot = curState.second.root;
		
		FA_DEBUG_AT(1,"Main processing product state (" << curState.first << ", " << curState.second << ") -> "
				<< curNewState.state << ", "
				<< "[" << thisRoot << "," << fwdRoot << "] -> " << curNewState.root);

		const std::shared_ptr<TreeAut> thisTA = thisFAE->getRoot(thisRoot);
		const std::shared_ptr<TreeAut> fwdTA = fwdFAE->getRoot(fwdRoot);
		assert((nullptr != thisTA) && (nullptr != fwdTA));

		const size_t& thisState = curState.first.state;
		const size_t& fwdState = curState.second.state;

		auto thisIt = thisTA->begin(thisState);
		auto thisEnd = thisTA->end(thisState);
		auto fwdEnd = fwdTA->end(fwdState);

		for (; thisIt != thisEnd; ++thisIt)
		{
			const auto& thisTrans = *thisIt;
			
			auto fwdIt = fwdTA->begin(fwdState);
			for (; fwdIt != fwdEnd; ++fwdIt)
			{
				const auto& fwdTrans = *fwdIt;

				// we handle data one level up
				assert(!TreeAut::GetSymbol(thisTrans)->isData()
                        && !TreeAut::GetSymbol(fwdTrans)->isData());
				assert(!thisTrans.GetChildren().empty()
                        && !fwdTrans.GetChildren().empty());

				// TODO: so far, we are not doing unfolding!
				if (TreeAut::GetSymbol(thisTrans) == TreeAut::GetSymbol(fwdTrans))
				{
					FA_DEBUG_AT(1,"Transition: " << TreeAut::GetSymbol(thisTrans) <<
                            ", thisState = " << thisState << ", srcState = " << fwdState);

					assert(thisTrans.GetChildrenSize() == fwdTrans.GetChildrenSize());
					const size_t& transArity = thisTrans.GetChildrenSize();

					if (TreeAut::GetSymbol(thisTrans)->isData())
					{	// data are processed one level up
						assert(TreeAut::GetSymbol(fwdTrans)->isData());
						continue;
					}

					// handle ordinary transitions
					std::vector<size_t> lhs;
					size_t i;
					for (i = 0; i < transArity; ++i)
					{	// for each pair of states that map to each other
						const Data* fwdData = nullptr, *thisData = nullptr;
						bool  fwdIsData =  fwdFAE->isData( fwdTrans.GetNthChildren(i),  fwdData);
						bool thisIsData = thisFAE->isData(thisTrans.GetNthChildren(i), thisData);

						if (!fwdIsData && !thisIsData)
						{	// ************* process internal states *************
							// This is the easiest case, when both states in the product are
							// internal. We do not create any new automaton.
							assert((nullptr == fwdData) && (nullptr == thisData));

							RootState rootState = engine.makeProductState(
								thisRoot, thisTrans.GetNthChildren(i),
								fwdRoot, fwdTrans.GetNthChildren(i));

							// check that we have not created a new automaton
							assert(rootState.root == curNewState.root);

							lhs.push_back(rootState.state);
						}
						else if (fwdIsData && thisIsData &&
							!fwdData->isRef() && !thisData->isRef())
						{ // ************* process real data states (leaves) *************
							// This is the second easiest case, when both states are real data
							// (i.e. no references). In this case, we are simply doing
							// intersection of the data.
							if (*thisData == *fwdData)
							{	// in case the data match
								lhs.push_back(fae->addData(
									*fae->getRoot(curNewState.root).get(), *thisData));
							}
							else
							{	// in case the data are different
								break;        // cut this branch of the intersection
							}
						}
						else if (fwdIsData && thisIsData &&
							fwdData->isRef() && thisData->isRef())
						{ // ************* process reference states (leaves) *************
							// This is another quite easy case, when both states are
							// references to another automata. In this case, we are jumping
							// from both automata into another product automaton.
							assert(thisData->d_ref.displ == fwdData->d_ref.displ);
							assert(0 == thisData->d_ref.displ);

							const size_t& thisNewRoot = thisData->d_ref.root;
							const size_t& fwdNewRoot  = fwdData->d_ref.root;

							const TreeAut* thisNewTA = thisFAE->getRoot(thisNewRoot).get();
							const TreeAut* fwdNewTA  = fwdFAE->getRoot(fwdNewRoot).get();
							assert((nullptr != thisNewTA) && (nullptr != fwdNewTA));

							RootState rootState = engine.makeProductState(
								thisNewRoot, thisNewTA->getFinalState(),
								fwdNewRoot, fwdNewTA->getFinalState());

							const size_t state = fae->addData(*fae->getRoot(curNewState.root).get(),
								Data::createRef(rootState.root));
							//fae->getRoot(curNewState.root)->addFinalState(state);
							FA_DEBUG_AT(1,"Adding data node r" << _MSB_GET(state) << " of root " << curNewState.root);
							lhs.push_back(state);
						}
						else if ((fwdIsData && !thisIsData && fwdData->isNull())
							|| (!fwdIsData && thisIsData && thisData->isNull())
							|| (fwdIsData && thisIsData && fwdData->isNull() && thisData->isRef())
							|| (fwdIsData && thisIsData && fwdData->isRef() && thisData->isNull()))
						{ // ************* process NULL pointers *************
							// This is the case when there is a NULL pointer and either an
							// internal state or a reference
							break;   // cut this branch of the intersection
						}
						else if ((fwdIsData && !thisIsData && fwdData->isRef())
							|| (!fwdIsData && thisIsData && thisData->isRef()))
						{ // ************* process jumps *************
							// This is the case when one FA jumps to another TA and the other does not
							FA_DEBUG_AT(1,"jump!");

							RootState rootState;
							if (fwdIsData)
							{
								assert(fwdIsData && !thisIsData && fwdData->isRef());

								const size_t& fwdNewRoot = fwdData->d_ref.root;
								const TreeAut* fwdNewTA  = fwdFAE->getRoot(fwdNewRoot).get();
								assert(nullptr != fwdNewTA);

								rootState = engine.makeProductState(
									thisRoot, thisTrans.GetNthChildren(i),
									fwdNewRoot, fwdNewTA->getFinalState());
							}
							else
							{
								assert(!fwdIsData && thisIsData && thisData->isRef());

								const size_t& thisNewRoot = thisData->d_ref.root;
								const TreeAut* thisNewTA = thisFAE->getRoot(thisNewRoot).get();
								assert(nullptr != thisNewTA);

								rootState = engine.makeProductState(
									thisNewRoot, thisNewTA->getFinalState(),
									fwdRoot, fwdTrans.GetNthChildren(i));
							}

							const size_t state = fae->addData(*fae->getRoot(curNewState.root).get(),
								Data::createRef(rootState.root));
							FA_DEBUG_AT(1,"Adding data node r" << _MSB_GET(state) << " of root " << curNewState.root);
							//fae->getRoot(curNewState.root)->addFinalState(state);
							lhs.push_back(state);
						}
						else if ((fwdIsData && fwdData->isUndef()) && thisIsData)
						{
							assert(fwdIsData && fwdData->isUndef());

							//FA_DEBUG_AT(1,"Substituting " << *srcData << " for " << newValue);

							lhs.push_back(fae->addData(
								*fae->getRoot(thisRoot).get(), *fwdData));
						}
						else
						{	// we should not get here
							assert(false);
						}
					}

					if (transArity == i)
					{	// in case we have not interrupted the search, add the transition
						std::ostringstream osLhs;
						for (auto it = lhs.cbegin(); it != lhs.cend(); ++it)
						{
							if (lhs.cbegin() != it)
								osLhs << ",";

							osLhs << FA::writeState(*it);
						}

						FA_DEBUG_AT(1,"TA " << curNewState.root << ": adding transition "
							<< FA::writeState(curNewState.state) << " -> "
							<< TreeAut::GetSymbol(thisTrans) << " (" << osLhs.str() << ")");

						fae->getRoot(curNewState.root)->addTransition(
							lhs, TreeAut::GetSymbol(thisTrans), curNewState.state);
					}
				}
			}
		}
	}

	FA_DEBUG_AT(1,"Result of intersection: " << *fae);

	// now, check whether there is some component with an empty language in the
	// result
	for (size_t i = 0; i < fae->getRootCount(); ++i)
	{
		TreeAut* ta = fae->allocTA();
		fae->getRoot(i)->uselessAndUnreachableFree(*ta);
		std::shared_ptr<TreeAut> pTa(ta);

		// check emptiness
		if (ta->getFinalStates().empty())
		{	// in case the language of an automaton is empty
			fae->clear();   // the language of the FA is empty
			FA_DEBUG_AT(1,"A tree of intersection is empty");
			return;
		}

		fae->setRoot(i, pTa);
	}

	// reorder the FAE to correspond to the original order
	std::vector<size_t> index = engine.getRootOrderIndexForRHS();
	assert(index.size() == fae->getValidRootCount());

	// set new root # for components that do not correspond to components in the
	// forward run


	std::ostringstream os;
	utils::printCont(os, index);
	FA_DEBUG_AT(1,"Index: " << os.str());

	std::vector<std::shared_ptr<TreeAut>> newRoots;
	size_t newRootsSize = std::max(fae->getRootCount(), fwdFAE->getRootCount());
	for (size_t i = 0; i < newRootsSize; ++i)
	{
		newRoots.push_back(std::shared_ptr<TreeAut>());
	}

	for (size_t i = 0; i < index.size(); ++i)
	{
		assert(index[i] < newRootsSize);
		newRoots[index[i]] = fae->getRoot(i);
	}

	// update representation
	fae->swapRoots(newRoots);

	FA_DEBUG_AT(1,"Before relabelling: " << *fae);

	for (size_t i = 0; i < index.size(); ++i)
	{
		fae->setRoot(index[i], std::shared_ptr<TreeAut>(
			fae->relabelReferences(fae->getRoot(index[i]).get(), index)
		));
	}

	FA_DEBUG_AT(1,"After shuffling: " << *fae);

	// FIXME: do we really need this?
//	fae->updateConnectionGraph();

	FA_DEBUG_AT(1,"Underapproximating intersection");
}


/*
void FAE::makeProduct(
	const FAE&                             lhs,
	const FAE&                             rhs,
	std::set<std::pair<size_t, size_t>>&   result)
{
	// engine that handles creation of new states etc.
	SubsRefEngine engine;

	assert(lhs.GetVarCount() == rhs.GetVarCount());
	for (size_t i = 0; i < lhs.GetVarCount(); ++i)
	{	// copy global variables
		const Data& lhsVar = lhs.GetVar(i);
		const Data& rhsVar = rhs.GetVar(i);

		if (lhsVar.isRef() && rhsVar.isRef())
		{	// for references
			assert(lhsVar.d_ref.displ == rhsVar.d_ref.displ);
			assert(0 == lhsVar.d_ref.displ);

			const TreeAut* lhsRoot = lhs.getRoot(lhsVar.d_ref.root).get();
			const TreeAut* rhsRoot = rhs.getRoot(rhsVar.d_ref.root).get();
			assert((nullptr != lhsRoot) && (nullptr != rhsRoot));

			engine.makeProductState(
				lhsVar.d_ref.root, lhsRoot->getFinalState(),
				rhsVar.d_ref.root, rhsRoot->getFinalState());
		}
	}

	while (!engine.wsEmpty())
	{
		const ProdState curState = engine.getNextState();

		const size_t& lhsRoot = curState.first.root;
		const size_t& rhsRoot = curState.second.root;

		const std::shared_ptr<TreeAut> lhsTA = lhs.getRoot(lhsRoot);
		const std::shared_ptr<TreeAut> rhsTA = rhs.getRoot(rhsRoot);
		assert((nullptr != lhsTA) && (nullptr != rhsTA));

		const size_t& lhsState = curState.first.state;
		const size_t& rhsState = curState.second.state;

		FA_DEBUG_AT(1,"Processing product state (" <<  curState.first << ", "
			<< curState.second << ")");

		auto lhsIt  = lhsTA->begin(lhsState);
		auto lhsEnd = lhsTA->end(lhsState);
		auto rhsEnd = rhsTA->end(rhsState);

		for (; lhsIt != lhsEnd; ++lhsIt)
		{
			auto rhsIt  = rhsTA->begin(rhsState);
			for (; rhsIt != rhsEnd; ++rhsIt)
			{
				const Transition& lhsTrans = *lhsIt;
				const Transition& rhsTrans = *rhsIt;

				// we handle data one level up
				assert(!TreeAut::GetSymbol(lhsTrans)->isData()
                        && !TreeAut::GetSymbol(rhsTrans)->isData());
				assert(!lhsTrans.GetChildren().empty() && !rhsTrans.GetChildren().empty());

				// TODO: so far, we are not doing unfolding!
				if (TreeAut::GetSymbol(lhsTrans) == TreeAut::GetSymbol(rhsTrans))
				{
					FA_DEBUG_AT(1,"Transition: " << TreeAut::GetSymbol(lhsTrans)
                            << ", thisState = " << lhsState << ", srcState = " << rhsState);

					assert(lhsTrans.GetChildrenSize() == rhsTrans.GetChildrenSize());

					if (TreeAut::GetSymbol(lhsTrans)->isData())
					{	// data are processed one level up
						assert(TreeAut::GetSymbol(rhsTrans)->isData());
						continue;
					}

					for (size_t i = 0; i < lhsTrans.GetChildrenSize(); ++i)
					{	// for each pair of states that map to each other
						const Data* lhsData = nullptr, *rhsData = nullptr;
						bool lhsIsData = lhs.isData(
                                lhsTrans.GetNthChildren(i), lhsData);
						bool rhsIsData = rhs.isData(
                                rhsTrans.GetNthChildren(i), rhsData);

						//  lhsIsData <-> nullptr != lrcData
						assert((!lhsIsData || (nullptr != lhsData))
							&& (lhsIsData || (nullptr == lhsData)));
						//  rhsIsData <-> nullptr != rhsData
						assert((!rhsIsData || (nullptr != rhsData))
							&& (rhsIsData || (nullptr == rhsData)));

						if (!lhsIsData && !rhsIsData)
						{	// ************* process internal states *************
							// This is the easiest case, when both states in the product are
							// internal.
							engine.makeProductState(
								lhsRoot, lhsTrans.GetNthChildren(i),
								rhsRoot, rhsTrans.GetNthChildren(i));
						}
						else if (lhsIsData && rhsIsData &&
							!lhsData->isRef() && !rhsData->isRef())
						{ // ************* process real data states (leaves) *************
							// This is the second easiest case, when both states are real data
							// (i.e. no references). In this case, we are simply doing
							// intersection of the data.
							if (*lhsData != *rhsData)
							{	// if data don't match
								break;        // cut this branch of the product
							}
						}
						else if (lhsIsData && rhsIsData &&
							lhsData->isRef() && rhsData->isRef())
						{ // ************* process reference states (leaves) *************
							// This is another quite easy case, when both states are
							// references to another automata. In this case, we are jumping
							// from both automata into another product automaton.
							assert(lhsData->d_ref.displ == rhsData->d_ref.displ);
							assert(0 == lhsData->d_ref.displ);

							FA_DEBUG_AT(1,"Two references");

							const size_t& lhsNewRoot = lhsData->d_ref.root;
							const size_t& rhsNewRoot = rhsData->d_ref.root;

							const TreeAut* lhsNewTA = lhs.getRoot(lhsNewRoot).get();
							const TreeAut* rhsNewTA = rhs.getRoot(rhsNewRoot).get();
							assert((nullptr != lhsNewTA) && (nullptr != rhsNewTA));

							engine.makeProductState(
								lhsNewRoot, lhsNewTA->getFinalState(),
								rhsNewRoot, rhsNewTA->getFinalState());
						}
						else if ((lhsIsData && !rhsIsData && lhsData->isNull())
							|| (!lhsIsData && rhsIsData && rhsData->isNull())
							|| (lhsIsData && rhsIsData && lhsData->isNull() && rhsData->isRef())
							|| (lhsIsData && rhsIsData && lhsData->isRef()  && rhsData->isNull()))
						{ // ************* process NULL pointers *************
							// This is the case when there is a NULL pointer and either an
							// internal state or a reference
							FA_DEBUG_AT(1,"NULL ptr!");

							break;   // cut this branch of the product
						}
						else if ((lhsIsData && !rhsIsData && lhsData->isRef())
							|| (!lhsIsData && rhsIsData && rhsData->isRef()))
						{ // ************* process jumps *************
							// This is the case when one FA jumps to another TA and the other does not
							FA_DEBUG_AT(1,"jump!");

							if (lhsIsData)
							{
								assert(lhsIsData && !rhsIsData && lhsData->isRef());

								const size_t& lhsNewRoot = lhsData->d_ref.root;
								const TreeAut* lhsNewTA  = lhs.getRoot(lhsNewRoot).get();
								assert(nullptr != lhsNewTA);

								engine.makeProductState(
									lhsNewRoot, lhsNewTA->getFinalState(),
									rhsRoot, rhsTrans.GetNthChildren(i));
							}
							else
							{
								assert(!lhsIsData && rhsIsData && rhsData->isRef());

								const size_t& rhsNewRoot = rhsData->d_ref.root;
								const TreeAut* rhsNewTA  = rhs.getRoot(rhsNewRoot).get();
								assert(nullptr != rhsNewTA);

								engine.makeProductState(
									lhsRoot, lhsTrans.GetNthChildren(i),
									rhsNewRoot, rhsNewTA->getFinalState());
							}
						}
						else
						{	// we should not get here
							assert(false);       // fail gracefully
						}
					}
				}
				else
				{
					FA_LOG("not-matching: " << *lhsIt << ", " << *rhsIt);
				}
			}
		}
	}

	for (const ProdState& prodState : engine.getProcessed())
	{
		result.insert(std::make_pair(prodState.first.state, prodState.second.state));
	}
}
 */
