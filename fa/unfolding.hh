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

#ifndef UNFOLDING_H
#define UNFOLDING_H

// Standard library headers
#include <vector>
#include <set>
#include <stdexcept>
#include <algorithm>

// Forester headers
#include "forestautext.hh"
#include "abstractbox.hh"
#include "utils.hh"

class Unfolding
{
	FAE& fae;

protected:

	/**
	 * @brief Merges two boxes to a new label
	 *
	 * Function takes two labels and joins them creating a new label.
	 * This is basically done during unfolding when a box is unfolded and
	 * the selectors to content of the box are added to
	 * a label of the wrapping automaton.
	 */
	void boxMerge(
		TreeAut&                       dst,
		const TreeAut&                 src,
		const TreeAut&                 boxRoot,
		const Box*                     box,
		const std::vector<size_t>&     rootIndex);

private:
    void getChildrenAndLabelFromBox(
		const Box*                         box,
        const TreeAut::Transition&         transition,
        std::vector<size_t>&               children,
        std::vector<const AbstractBox*>&   label);

	void initRootRefIndex(
		std::vector<size_t>&                   index,
		const Box*                             box,
		const TreeAut::Transition&             t);

	void substituteOutputPorts(
		const std::vector<size_t>&    index,
		const size_t                  root,
		const Box*                    box);

	void substituteInputPorts(
		const std::vector<size_t>&    index,
		const Box*                    box);

public:

	typedef std::unordered_map<const Box*, size_t>		BoxSignature;
	typedef std::unordered_map<size_t, BoxSignature>	StateToBoxSignatureMap;

private:

	static bool updateBoxSignature(
		StateToBoxSignatureMap&	stateMap,
		size_t			state,
		const BoxSignature&	boxSignature);

	static void joinBoxSignature(
		BoxSignature&				dst,
		const std::pair<const Box*, size_t>&	boxOriginPair);

	static void joinBoxSignature(BoxSignature& dst, const BoxSignature& src);

	static bool processLhs(
		BoxSignature&			result,
		const std::vector<size_t>&	lhs,
		const StateToBoxSignatureMap&	stateMap);

public:

	static void computeBoxSignatureMap(
		StateToBoxSignatureMap&	stateMap,
		const TreeAut&		ta);

public:
	Unfolding(FAE& fae) : fae(fae) {}

	void unfoldBox(const size_t root, const Box* box);
	void unfoldBoxes(const size_t root, const std::set<const Box*>& boxes);

	void unfoldSingletons(const size_t root);
};

#endif
