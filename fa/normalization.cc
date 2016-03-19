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

// Code Listener headers
#include <cl/storage.hh>

// Forester headers
#include "abstractinstruction.hh"
#include "error_messages.hh"
#include "garbage_checker.hh"
#include "normalization.hh"
#include "programerror.hh"
#include "regdef.hh"
#include "symstate.hh"
#include "virtualmachine.hh"

namespace {
	
	void substituteRootRefByTA(
		TreeAut&                        resTA,
		const TreeAut::Transition&      trans,
		const size_t                    rootRefState,
		const std::vector<size_t>&      joinStates)
	{
		std::vector<size_t> transChildren = trans.GetChildren();
		std::vector<size_t>::iterator rootRefIterator = 
			std::find(transChildren.begin(), transChildren.end(), rootRefState);
		assert(rootRefIterator != transChildren.end());

		for (const auto& k : joinStates)
		{
			*rootRefIterator = k;
			resTA.addTransition(
				transChildren,
				TreeAut::GetSymbol(trans),
				trans.GetParent());
		}
	}

	// StateTuple has to be iterable and has to have defined equivalence relation
	template<class StateTuple>
	static bool containsRootRef(
		const StateTuple& stateTuple,
		const size_t rootRefState)
	{
		return std::find(stateTuple.begin(), stateTuple.end(), rootRefState) != stateTuple.end();
	}
}

TreeAut* Normalization::mergeRoot(
	TreeAut&                          dst,
	const size_t                      rootRef,
	const TreeAut&                    src,
	std::vector<size_t>&              joinStates)
{
	// Assertions
	assert(rootRef < this->fae.getRootCount());

	std::unordered_map<size_t, size_t> joinStatesMap;
	for (const size_t finState : src.getFinalStates())
	{
		joinStates.push_back(this->fae.nextState());
		joinStatesMap.insert(std::make_pair(finState, this->fae.freshState()));
	}

	TreeAut* resTA = this->fae.allocTA();
	resTA->addFinalStates(dst.getFinalStates());
	size_t rootRefState = _MSB_ADD(this->fae.boxMan->getDataId(Data::createRef(rootRef)));
	
	bool hit = false;
	for (const auto& trans : dst)
	{
		if (containsRootRef(trans.GetChildren(), rootRefState))
		{
			substituteRootRefByTA(*resTA, trans, rootRefState, joinStates);
			hit = true;
		}
		else
		{
			resTA->addTransition(trans);
		}
	}
	if (!hit) {assert(false);}
	// avoid screwing up things
	src.unfoldAtRoot(*resTA, joinStatesMap, false);

	// postconditions
	
	// the new TA should have the same final states as dst TA
	// so at least number of the states is checked. 
	assert(resTA->getFinalStates().size() == dst.getFinalStates().size());
	return resTA;
}


void Normalization::traverse(
	std::vector<bool>&                visited,
	std::vector<size_t>&              order,
	std::vector<bool>&                marked) const
{
	visited = std::vector<bool>(this->fae.getRootCount(), false);
	marked = std::vector<bool>(this->fae.getRootCount(), false);

	order.clear();

	for (const Data& var : this->fae.GetVariables())
	{	// start from program variables and perform DFT
		// skip everything what is not a root reference
		if (!var.isRef())
			continue;

		size_t root = var.d_ref.root;

		// mark rootpoint pointed by a variable
		marked[root] = true;

		// check whether we traversed this one before
		if (visited[root])
			continue;

		this->fae.connectionGraph.visit(root, visited, order, marked);
	}
}


void Normalization::normalizeRoot(
	std::vector<bool>&                normalized,
	const size_t                      root,
	const std::vector<bool>&          marked)
{
	if (normalized[root])
		return;

	normalized[root] = true;

	// we need a copy here!
	ConnectionGraph::CutpointSignature signature =
		this->fae.connectionGraph.data[root].signature;

	for (auto& cutpoint : signature)
	{
		this->normalizeRoot(normalized, cutpoint.root, marked);

		if (marked[cutpoint.root])
			continue;

		assert(root != cutpoint.root);

		std::vector<size_t> refStates;

		TreeAut* ta = this->mergeRoot(
			*this->fae.getRoot(root),
			cutpoint.root,
			*this->fae.getRoot(cutpoint.root),
			refStates
		);

		this->fae.setRoot(root, std::shared_ptr<TreeAut>(ta));
		this->fae.setRoot(cutpoint.root, nullptr);

		this->fae.connectionGraph.mergeCutpoint(root, cutpoint.root);
	}
}


bool Normalization::selfReachable(
	size_t                            root,
	size_t                            self,
	const std::vector<bool>&          marked)
{
	for (auto& cutpoint : this->fae.connectionGraph.data[root].signature)
	{
		if (cutpoint.root == self)
			return true;

		if (marked[cutpoint.root])
			continue;

		if (this->selfReachable(cutpoint.root, self, marked))
			return true;
	}

	return false;
}


void Normalization::scan(
	std::vector<bool>&                marked,
	std::vector<size_t>&              order,
	const std::set<size_t>&           forbidden,
	bool                              extended)
{
	assert(this->fae.connectionGraph.isValid());

	std::vector<bool> visited(this->fae.getRootCount(), false);

	marked = std::vector<bool>(this->fae.getRootCount(), false);

	order.clear();

	// compute canonical root ordering
	this->traverse(visited, order, marked);

	GarbageChecker::nontraverseCheckAndRemoveGarbage(this->fae, this->state_, visited);

	if (!extended)
	{
		for (auto& x : forbidden)
		{
			marked[x] = true;
		}

		return;
	}

	for (auto& x : forbidden)
	{
		marked[x] = true;

		for (auto& cutpoint : this->fae.connectionGraph.data[x].signature)
		{
			if ((cutpoint.root != x) && !this->selfReachable(cutpoint.root, x, marked))
				continue;

			marked[cutpoint.root] = true;

			break;
		}
	}
}

namespace {
	static bool isInCanonicalForm(
		const std::vector<bool>&          marked,
		const std::vector<size_t>&        order)
	{
		size_t i;

		for (i = 0; i < order.size(); ++i)
		{
			if (!marked[i] || (order[i] != i))
				return false;
		}

		assert(i == order.size());
		return true;
	}
}

bool Normalization::normalizeInternal(
	const std::vector<bool>&          marked,
	const std::vector<size_t>&        order)
{
	if (isInCanonicalForm(marked, order))
	{	// in case the FA is in the canonical form
		this->fae.resizeRoots(order.size());
		this->fae.connectionGraph.data.resize(order.size());
		return false;
	}

	// reindex roots
	std::vector<size_t> index(this->fae.getRootCount(), static_cast<size_t>(-1));
	std::vector<bool> normalized(this->fae.getRootCount(), false);
	std::vector<std::shared_ptr<TreeAut>> newRoots;
	size_t offset = 0;
	bool merged = false;

	for (auto& i : order)
	{	// push tree automata into a new tuple in the right order
		this->normalizeRoot(normalized, i, marked);

		if (!marked[i])
		{	// if a root was merged, do not put it in the new tuple!
			merged = true;

			continue;
		}

		newRoots.push_back(this->fae.getRoot(i));

		index[i] = offset++;
	}

	// update representation
	this->fae.swapRoots(newRoots);

	for (size_t i = 0; i < this->fae.getRootCount(); ++i)
	{
		this->fae.setRoot(i, std::shared_ptr<TreeAut>(
			this->fae.relabelReferences(this->fae.getRoot(i).get(), index)
		));
	}

	this->fae.connectionGraph.finishNormalization(this->fae.getRootCount(), index);

	// update variables
	this->fae.UpdateVarsRootRefs(index);

	return merged;
}

bool Normalization::normalize(
		FAE&                              fae,
		const SymState*                   state,
		const std::set<size_t>&           forbidden,
		bool                              extended)
{
	Normalization norm(fae, state);

	std::vector<size_t> order;
	std::vector<bool> marked;

	norm.scan(marked, order, forbidden, extended);

	bool result = norm.normalizeInternal(marked, order);

	FA_DEBUG_AT(3, "after normalization: " << std::endl << fae);

	return result;
}

bool Normalization::normalizeWithoutMerging(
		FAE&                              fae,
		const SymState*                   state,
		const std::set<size_t>&           forbidden,
		bool                              extended)
{
	Normalization norm(fae, state);

	std::vector<size_t> order;
	std::vector<bool> marked;

	norm.scan(marked, order, forbidden, extended);

	// normalize without merging (we say that all components are referred more
	// than once), i.e. only reorder
	std::fill(marked.begin(), marked.end(), true);
	
	bool result = norm.normalizeInternal(marked, order);

	FA_DEBUG_AT(3, "after normalization: " << std::endl << fae);

	return result;
}

std::set<size_t> Normalization::computeForbiddenSet(FAE& fae)
{
	// Assertions
	assert(fae.getRootCount() == fae.connectionGraph.data.size());
	assert(fae.getRootCount() >= FIXED_REG_COUNT);

	std::set<size_t> forbidden;

	VirtualMachine vm(fae);

	for (size_t i = 0; i < FIXED_REG_COUNT; ++i)
	{
		assert(fae.getRoot(vm.varGet(i).d_ref.root));
		forbidden.insert(vm.varGet(i).d_ref.root);
	}

	for (size_t i = 0; i < FIXED_REG_COUNT; ++i)
	{
		vm.getNearbyReferences(vm.varGet(i).d_ref.root, forbidden);
	}

	return forbidden;
}
