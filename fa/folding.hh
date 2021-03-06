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

#ifndef FOLDING_H
#define FOLDING_H

// Standard library headers
#include <vector>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

// Forester headers
#include "abstractbox.hh"
#include "boxman.hh"
#include "config.h"
#include "connection_graph.hh"
#include "forestautext.hh"
#include "restart_request.hh"
#include "unfolding.hh"
#include "streams.hh"

class Folding {

protected:

	static void copyBox(std::vector<size_t>& lhs, std::vector<const AbstractBox*>& label,
		const AbstractBox* box, const std::vector<size_t>& srcLhs, const size_t& srcOffset) {

		for (size_t i = 0; i < box->getArity(); ++i)
			lhs.push_back(srcLhs[srcOffset + i]);

		label.push_back(box);

	}

	static const ConnectionGraph::CutpointSignature& getSignature(size_t state,
		const ConnectionGraph::StateToCutpointSignatureMap& signatures) {

		auto iter = signatures.find(state);

		if (iter == signatures.end())
			FA_DEBUG_AT(3, state);

		assert(iter != signatures.end());

		return iter->second;

	}

	static bool isSignaturesCompatible(const ConnectionGraph::CutpointSignature& s1,
		const ConnectionGraph::CutpointSignature& s2) {

		if (s1.size() != s2.size())
			return false;

		for (size_t i = 0; i < s1.size(); ++i) {

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

	const ConnectionGraph::StateToCutpointSignatureMap& getSignatures(size_t root) {

		assert(root < this->signatureMap.size());

		if (!this->signatureMap[root].first) {

			ConnectionGraph::computeSignatures(
				this->signatureMap[root].second, *this->fae.roots[root]
			);

			this->signatureMap[root].first = true;

		}

		assert(this->signatureMap[root].first);

		return this->signatureMap[root].second;

	}

	void invalidateSignatures(size_t root) {

		assert(root < this->signatureMap.size());

		this->signatureMap[root].first = false;

	}

	bool componentCut(TreeAut& res, TreeAut& complement,
		ConnectionGraph::CutpointSignature& complementSignature, size_t root, size_t state,
		size_t target
	) {

		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		const TreeAut& src = *this->fae.roots[root];

		res.addFinalStates(src.getFinalStates());

		complement.addFinalState(state);

		auto& signatures = this->getSignatures(root);

		std::unordered_set<const AbstractBox*> boxes;

		// enumerate which boxes to fold
		for (auto i = src.begin(state); i != src.end(state, i); ++i) {

			size_t lhsOffset = 0;

			const TT<label_type>& t = *i;

			for (auto& box : t.label()->getNode()) {

				// look for target cutpoint
				for (size_t j = 0; j < box->getArity(); ++j) {

					assert(lhsOffset + j < t.lhs().size());

					if (ConnectionGraph::containsCutpoint(
						Folding::getSignature(t.lhs()[lhsOffset + j], signatures), target)
					)
					{
						boxes.insert(box);

						break;
					}

				}

				lhsOffset += box->getArity();

			}

		}

		ConnectionGraph::CutpointSignature tmp;

		for (auto i = src.begin(); i != src.end(); ++i) {

			const TT<label_type>& t = *i;

			if (t.rhs() != state) {

				res.addTransition(t);
				complement.addTransition(t);

				continue;

			}

			std::vector<size_t> lhs, cLhs;
			std::vector<const AbstractBox*> label, cLabel;

			size_t lhsOffset = 0;

			tmp.clear();

			// split transition
			for (auto& box : t.label()->getNode()) {

				if (boxes.count(box) == 0) {

					Folding::copyBox(lhs, label, box, t.lhs(), lhsOffset);

				} else {

					for (size_t j = 0; j < box->getArity(); ++j) {

						assert(lhsOffset + j < t.lhs().size());

						ConnectionGraph::processStateSignature(
							tmp,
							static_cast<const StructuralBox*>(box),
							j,
							t.lhs()[lhsOffset + j],
							Folding::getSignature(t.lhs()[lhsOffset + j], signatures)
						);

					}

					Folding::copyBox(cLhs, cLabel, box, t.lhs(), lhsOffset);

				}

				lhsOffset += box->getArity();

			}

			ConnectionGraph::normalizeSignature(tmp);

			assert(tmp.size());

			// did we arrive here for the first time?
			if (complementSignature.empty())
				complementSignature = tmp;

			if (!Folding::isSignaturesCompatible(complementSignature, tmp))
			{	// seems that cut is not possible here
				FA_DEBUG_AT(2, "cannot cut: " << complementSignature << " != " << tmp);

				return false;
			}

			for (size_t i = 0; i < tmp.size(); ++i) {

				complementSignature[i].refCount = std::max(complementSignature[i].refCount, tmp[i].refCount);

				complementSignature[i].fwdSelectors.insert(
					tmp[i].fwdSelectors.begin(), tmp[i].fwdSelectors.end()
				);

			}

			assert(label.size());
			FAE::reorderBoxes(label, lhs);
			res.addTransition(lhs, this->fae.boxMan->lookupLabel(label), state);

			assert(cLabel.size());
			FAE::reorderBoxes(cLabel, cLhs);
			complement.addTransition(cLhs, this->fae.boxMan->lookupLabel(cLabel), state);

		}

		return true;

	}

	bool separateCutpoint(
		std::pair<std::shared_ptr<TreeAut>, std::shared_ptr<TreeAut>>& result,
		ConnectionGraph::CutpointSignature& boxSignature, size_t root, size_t state,
		size_t cutpoint)
	{
		auto ta = std::shared_ptr<TreeAut>(this->fae.allocTA());
		auto tmp = std::shared_ptr<TreeAut>(this->fae.allocTA());

		if (!this->componentCut(*ta, *tmp, boxSignature, root, state, cutpoint))
			return false;

		auto tmp2 = std::shared_ptr<TreeAut>(this->fae.allocTA());

		tmp->unreachableFree(*tmp2);

		result = std::make_pair(ta, tmp2);

		return true;
	}

	std::shared_ptr<TreeAut> relabelReferences(const TreeAut& ta,
		std::vector<size_t>& index) {

		auto tmp = std::shared_ptr<TreeAut>(this->fae.allocTA());

		this->fae.relabelReferences(*tmp, ta, index);

		return tmp;

	}

	/**
	 * @brief  @todo
	 */
	std::shared_ptr<TreeAut> joinBox(const TreeAut& src, size_t state, size_t root,
		const Box* box, const ConnectionGraph::CutpointSignature& signature) {

		auto ta = std::shared_ptr<TreeAut>(this->fae.allocTA());

		ta->addFinalStates(src.getFinalStates());

		for (auto i = src.begin(); i != src.end(); ++i) {

			if (i->rhs() != state) {

				ta->addTransition(*i);

				continue;

			}

			std::vector<const AbstractBox*> label(i->label()->getNode());
			std::vector<size_t> lhs(i->lhs());

			label.push_back(box);

			for (auto& cutpoint : signature) {

				if ((cutpoint.root == root) && (src.getFinalStates().count(state)))
					continue;

				lhs.push_back(
					this->fae.addData(*ta, Data::createRef(cutpoint.root))
				);

			}

			FA::reorderBoxes(label, lhs);

			ta->addTransition(lhs, this->fae.boxMan->lookupLabel(label), state);

		}

		return ta;

	}

	static void updateSelectorMap(std::unordered_map<size_t, size_t>& m, size_t selector,
		const ConnectionGraph::CutpointSignature& signature) {

		for (auto& cutpoint : signature) {

			auto p = m.insert(std::make_pair(cutpoint.root, selector));

			if (!p.second && p.first->second > selector)
				p.first->second = selector;

		}

	}

	// compute cutpoint-to-selector mapping, i.e. tell which selector one needs to take
	// in order to reach a given cutpoint
	static void computeSelectorMap(std::unordered_map<size_t, size_t>& selectorMap,
		const TT<label_type>& t, const ConnectionGraph::StateToCutpointSignatureMap& stateMap) {

		size_t lhsOffset = 0;

		for (auto& absBox : t.label()->getNode()) {

			switch (absBox->getType()) {

				case box_type_e::bSel: {

					auto iter = stateMap.find(t.lhs()[lhsOffset]);

					assert(iter != stateMap.end());

					Folding::updateSelectorMap(
						selectorMap, (static_cast<const SelBox*>(absBox))->getData().offset, iter->second
					);

					break;
				}

				case box_type_e::bBox: {

					const Box* box = static_cast<const Box*>(absBox);

					for (size_t i = 0; i < box->getArity(); ++i) {

						auto iter = stateMap.find(t.lhs()[lhsOffset + i]);

						assert(iter != stateMap.end());

						Folding::updateSelectorMap(
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

	bool checkSelectorMap(const std::unordered_map<size_t, size_t>& selectorMap,
		size_t root, size_t state) {

		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		auto& signatures = this->getSignatures(root);

		auto& ta = *this->fae.roots[root];

		for (TreeAut::iterator i = ta.begin(state); i != ta.end(state, i); ++i) {

			std::unordered_map<size_t, size_t> m;

			Folding::computeSelectorMap(m, *i, signatures);

			if (m != selectorMap)
				return false;

		}

		return true;

	}

	bool computeSelectorMap(std::unordered_map<size_t, size_t>& selectorMap,
		size_t root, size_t state) {

		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		auto& ta = *this->fae.roots[root];

		assert(ta.begin(state) != ta.end(state));
/*
		for (TreeAut::iterator i = ta.accBegin(); i != ta.accEnd(i); ++i)
			Folding::computeSelectorMap(selectorMap, *i, stateMap);
*/
		auto& signatures = this->getSignatures(root);

		Folding::computeSelectorMap(selectorMap, *ta.begin(state), signatures);

		return this->checkSelectorMap(selectorMap, root, state);

	}

	static size_t extractSelector(const std::unordered_map<size_t, size_t>& selectorMap,
		size_t target) {

		auto iter = selectorMap.find(target);

		assert(iter != selectorMap.end());

		return iter->second;

	}

	// transform
	static void extractInputMap(std::vector<size_t>& inputMap,
		const std::unordered_map<size_t, size_t>& selectorMap, size_t root,
		const std::vector<size_t>& index) {

		assert(index[root] == 0);

		inputMap.resize(selectorMap.size());

		size_t count = 0;

		for (auto& cutpointSelectorPair : selectorMap) {

			if (cutpointSelectorPair.first == root) {

				// reference to root does not appear in the box interface
				continue;

			}

			assert(cutpointSelectorPair.first < index.size());

			if (index[cutpointSelectorPair.first] == static_cast<size_t>(-1))
				continue;

			assert(index[cutpointSelectorPair.first] >= 1);
			assert(index[cutpointSelectorPair.first] < inputMap.size() + 1);

			inputMap[index[cutpointSelectorPair.first] - 1] = cutpointSelectorPair.second;

			++count;

		}

		inputMap.resize(count);

	}
/*
	static bool checkSingular(const TreeAut& ta, bool result,
		const ConnectionGraph::StateToCutpointSignatureMap& stateMap) {

		for (auto& state : ta.getFinalStates()) {

			auto iter = stateMap.find(state);

			assert(iter != stateMap.end());

			if (iter->second.empty() != result)
				return false;

		}

		return true;

	}

	static bool isSingular(const TreeAut& ta) {

		ConnectionGraph::StateToCutpointSignatureMap stateMap;

		ConnectionGraph::computeSignatures(stateMap, ta);

		assert(ta.getFinalStates().size());

		auto iter = stateMap.find(*ta.getFinalStates().begin());

		assert(iter != stateMap.end());

		bool result = iter->second.empty();

		assert(Folding::checkSingular(ta, result, stateMap));

		return result;

	}
*/
	const Box* getBox(const Box& box, bool conditional) {

		return (conditional)?(this->boxMan.lookupBox(box)):(this->boxMan.getBox(box));

	}

	const Box* makeType1Box(size_t root, size_t state, size_t aux, const std::set<size_t>& forbidden,
		bool conditional = true, bool test = false) {

		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		std::vector<size_t> index(this->fae.roots.size(), static_cast<size_t>(-1)), inputMap;
		std::unordered_map<size_t, size_t> selectorMap;
		ConnectionGraph::CutpointSignature outputSignature;

		size_t start = 0;

		std::pair<std::shared_ptr<TreeAut>, std::shared_ptr<TreeAut>> p;

		if (!this->separateCutpoint(p, outputSignature, root, state, aux))
			return nullptr;

		index[root] = start++;

		for (auto& cutpoint : outputSignature) {

			if (forbidden.count(cutpoint.root))
				return nullptr;

			assert(cutpoint.root < index.size());

			if (cutpoint.root != root)
				index[cutpoint.root] = start++;

		}

		if (!Folding::computeSelectorMap(selectorMap, root, state))
			return nullptr;

		Folding::extractInputMap(inputMap, selectorMap, root, index);

		auto box = std::unique_ptr<Box>(
			this->boxMan.createType1Box(
				root,
				this->relabelReferences(*p.second, index),
				outputSignature,
				inputMap,
				index
			)
		);

		auto boxPtr = this->getBox(*box, conditional);

		if (test)
			return boxPtr;

		if (!boxPtr)
			return nullptr;

		FA_DEBUG_AT(2, *static_cast<const AbstractBox*>(boxPtr) << " found");

		this->fae.roots[root] = this->joinBox(*p.first, state, root, boxPtr, outputSignature);
		this->fae.connectionGraph.invalidate(root);

		this->invalidateSignatures(root);

		return boxPtr;

	}

	const Box* makeType2Box(size_t root, size_t aux, const std::set<size_t>& forbidden,
		bool conditional = true, bool test = false) {

		assert(root < this->fae.roots.size());
		assert(aux < this->fae.roots.size());
		assert(this->fae.roots[root]);
		assert(this->fae.roots[aux]);

		size_t finalState = this->fae.roots[root]->getFinalState();

		std::vector<size_t> index(this->fae.roots.size(), static_cast<size_t>(-1)), index2, inputMap;
		std::vector<bool> rootMask(this->fae.roots.size(), false);
		std::unordered_map<size_t, size_t> selectorMap;
		ConnectionGraph::CutpointSignature outputSignature, inputSignature, tmpSignature;

		size_t start = 0;

		std::pair<std::shared_ptr<TreeAut>, std::shared_ptr<TreeAut>> p;

		if (!this->separateCutpoint(p, outputSignature, root, finalState, aux))
			return nullptr;

		index[root] = start++;

		for (auto& cutpoint : outputSignature) {
/*
			// we assume type 1 box is not present
			assert(cutpoint.root != root);
*/
			if (cutpoint.root == root)
				return nullptr;

			if (forbidden.count(cutpoint.root))
				return nullptr;

			assert(cutpoint.root < index.size());

			if (cutpoint.root != root)
				index[cutpoint.root] = start++;

		}

		if (!Folding::computeSelectorMap(selectorMap, root, finalState))
			return nullptr;

		Folding::extractInputMap(inputMap, selectorMap, root, index);

		std::pair<std::shared_ptr<TreeAut>, std::shared_ptr<TreeAut>> auxP;

		if (!this->separateCutpoint(auxP, inputSignature, aux, this->fae.roots[aux]->getFinalState(), root))
			return nullptr;
/*
		if (Folding::isSingular(*auxP.first))
			return false;
*/
		index2 = index;

		for (auto& cutpoint : inputSignature) {

			if (cutpoint.refCount > 1)
				return nullptr;

			if (forbidden.count(cutpoint.root))
				return nullptr;

			assert(cutpoint.root < index.size());

			if (index[cutpoint.root] == static_cast<size_t>(-1)) {

				assert(index2[cutpoint.root] == static_cast<size_t>(-1));

				index2[cutpoint.root] = start++;

				tmpSignature.push_back(cutpoint);

				inputMap.push_back(static_cast<size_t>(-1));

			}

		}

		selectorMap.clear();

		if (!Folding::computeSelectorMap(selectorMap, aux, this->fae.roots[aux]->getFinalState()))
			assert(false);

		size_t selector = Folding::extractSelector(selectorMap, root);

		auto box = std::unique_ptr<Box>(
			this->boxMan.createType2Box(
				root,
				this->relabelReferences(*p.second, index),
				outputSignature,
				inputMap,
				aux,
				this->relabelReferences(*auxP.second, index2),
				inputSignature,
				selector,
				index
			)
		);

		auto boxPtr = this->getBox(*box, conditional);

		if (test)
			return boxPtr;

		if (!boxPtr)
			return nullptr;

		FA_DEBUG_AT(2, *static_cast<const AbstractBox*>(boxPtr) << " found");

		for (auto& cutpoint : tmpSignature)
			outputSignature.push_back(cutpoint);

		this->fae.roots[root] = this->joinBox(*p.first, finalState, root, boxPtr, outputSignature);
		this->fae.connectionGraph.invalidate(root);

		this->invalidateSignatures(root);

		this->fae.roots[aux] = auxP.first;
		this->fae.connectionGraph.invalidate(aux);

		this->invalidateSignatures(aux);

		return boxPtr;

	}

public:

	bool discover1(size_t root, const std::set<size_t>& forbidden, bool conditional) {

		assert(this->fae.roots.size() == this->fae.connectionGraph.data.size());
		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		if (forbidden.count(root))
			return false;

		this->fae.updateConnectionGraph();

		if (!ConnectionGraph::containsCutpoint(this->fae.connectionGraph.data[root].signature, root))
			return false;

#if FA_TYPE_1_UNFOLD_HEURISTICS
		if (this->fae.roots[root]->getAcceptingTransitionCount() == 1)
		{
			const TT<label_type>& t = this->fae.roots[root]->getAcceptingTransition();

			std::set<const Box*> boxes;

			for (auto aBox : t.label()->getNode())
			{
				if ((aBox->getArity() == 0) && aBox->isBox() && static_cast<const Box*>(aBox)->hasSelfReference())
					boxes.insert(static_cast<const Box*>(aBox));
			}

			if (boxes.size())
			{
				FA_DEBUG_AT(3, "unfolding type 1 box at root " << root);

				Unfolding(this->fae).unfoldBoxes(root, boxes);

				this->signatureMap[root].first = false;

				this->fae.updateConnectionGraph();

				FA_DEBUG_AT(3, "after unfolding: " << std::endl << this->fae);
			}
		}
#endif

		bool found = false, hit;

		do
		{
			FA_DEBUG_AT(3, "type 1 cutpoint detected at root " << root);

			// save state offset
			this->fae.pushStateOffset();

			hit = (this->makeType1Box(
				root, this->fae.roots[root]->getFinalState(), root, forbidden, conditional
				) != nullptr);

			if (hit)
			{
				found = true;

				this->signatureMap[root].first = false;

				this->fae.updateConnectionGraph();
			}

		} while (hit && ConnectionGraph::containsCutpoint(this->fae.connectionGraph.data[root].signature, root));

		this->fae.popStateOffset();

		return found;

	}

	bool discover2(size_t root, const std::set<size_t>& forbidden, bool conditional) {

		assert(this->fae.roots.size() == this->fae.connectionGraph.data.size());
		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		if (forbidden.count(root))
			return nullptr;

		bool found = false;
dis2_start:
		// save state offset
		this->fae.pushStateOffset();

		this->fae.updateConnectionGraph();

		for (auto& cutpoint : this->fae.connectionGraph.data[root].signature) {

			if (cutpoint.refCount < 2)
				continue;

			auto& signatures = this->getSignatures(root);

			for (auto& stateSignaturePair : signatures) {

				for (auto& tmp : stateSignaturePair.second) {

					if ((tmp.refCount < 2) || tmp.refInherited || (tmp.root == cutpoint.root))
						continue;

					FA_DEBUG_AT(3, "type 2 cutpoint detected inside component " << root << " at state q" << stateSignaturePair.first);

					auto boxPtr = this->makeType1Box(
						root, stateSignaturePair.first, cutpoint.root, forbidden, conditional
					);

					if (boxPtr) {

						found = true;

						this->signatureMap[root].first = false;

						goto dis2_start;

					}

					this->fae.popStateOffset();

				}

			}

		}

		return found;

	}

	bool discover3(size_t root, const std::set<size_t>& forbidden, bool conditional) {

		assert(this->fae.roots.size() == this->fae.connectionGraph.data.size());
		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		if (forbidden.count(root))
			return nullptr;

		bool found = false;
dis3_start:
		// save state offset
		this->fae.pushStateOffset();

		this->fae.updateConnectionGraph();

		for (auto& cutpoint : this->fae.connectionGraph.data[root].signature) {

			if (forbidden.count(cutpoint.root)/* || cutpoint.joint*/)
				continue;

			size_t selectorToRoot = ConnectionGraph::getSelectorToTarget(
				this->fae.connectionGraph.data[cutpoint.root].signature, root
			);

			if (selectorToRoot == static_cast<size_t>(-1))
				continue;
/*
			if (selectorToRoot == cutpoint.forwardSelector)
				continue;
*/
			assert(!cutpoint.fwdSelectors.empty());

			if (/*(selectorToRoot < *cutpoint.fwdSelectors.begin()) ||*/
				this->makeType2Box(cutpoint.root, root, forbidden, true, true))
					continue;

			FA_DEBUG_AT(3, "type 3 cutpoint detected at roots " << root << " and " << cutpoint.root);

			auto boxPtr = this->makeType2Box(root, cutpoint.root, forbidden, conditional);

			if (boxPtr) {

				found = true;

				this->signatureMap[root].first = false;
				this->signatureMap[cutpoint.root].first = false;

				goto dis3_start;

			}

			this->fae.popStateOffset();

		}

		return found;

	}
/*
	bool discover(size_t root, const std::set<size_t>& forbidden) {

		assert(this->fae.connectionGraph.isValid());
		assert(this->fae.roots.size() == this->fae.connectionGraph.data.size());
		assert(root < this->fae.roots.size());
		assert(this->fae.roots[root]);

		if (forbidden.count(root))
			return false;

		const Box* boxPtr;

		FA_DEBUG_AT(3, "analysing: " << this->fae);

		// save state offset
		this->fae.pushStateOffset();

		bool found = false;

		size_t selectorToRoot;
start:
		this->fae.updateConnectionGraph();

		for (auto& cutpoint : this->fae.connectionGraph.data[root].signature) {

			if (cutpoint.root == root) {

				FA_DEBUG_AT(3, "type 1 cutpoint detected at root " << root);

				boxPtr = this->makeType1Box(
					root, this->fae.roots[root]->getFinalState(), root, forbidden, true
				);

				if (boxPtr)
					goto box_found;

				this->fae.popStateOffset();

				continue;

			}

			if (cutpoint.joint) {

				auto& signatures = this->getSignatures(root);

				for (auto& stateSignaturePair : signatures) {

					for (auto& tmp : stateSignaturePair.second) {

						if (!tmp.joint || tmp.joinInherited || (tmp.root != cutpoint.root))
							continue;

						FA_DEBUG_AT(3, "type 2 cutpoint detected inside component " << root << " at state q" << stateSignaturePair.first);

						boxPtr = this->makeType1Box(
							root, stateSignaturePair.first, cutpoint.root, forbidden, true
						);

						if (boxPtr)
							goto box_found;

						this->fae.popStateOffset();

					}

				}

				continue;

			}

			if (forbidden.count(cutpoint.root))
				continue;

			selectorToRoot = ConnectionGraph::getSelectorToTarget(
				this->fae.connectionGraph.data[cutpoint.root].signature, root
			);

			if (selectorToRoot == (size_t)(-1))
				continue;
*//*
			if (selectorToRoot == cutpoint.forwardSelector)
				continue;
*//*
			assert(!cutpoint.fwdSelectors.empty());

			FA_DEBUG_AT(3, "type 3 cutpoint detected at roots " << root << " and " << cutpoint.root);

			boxPtr = this->makeType2Box(root, cutpoint.root, forbidden, true);

			if (boxPtr)
				goto box_found;

			this->fae.popStateOffset();

			continue;
box_found:
			FA_DEBUG_AT(3, (AbstractBox*)boxPtr << " found");

			found = true;

			if (!this->fae.connectionGraph.data[root].valid)
				goto start;

		}

		return found;

	}
*/
public:

	Folding(FAE& fae, BoxMan& boxMan) : fae(fae), boxMan(boxMan), signatureMap(fae.roots.size()) {}

private:

	FAE& fae;
	BoxMan& boxMan;

	std::vector<std::pair<bool, ConnectionGraph::StateToCutpointSignatureMap>> signatureMap;

};

#endif
