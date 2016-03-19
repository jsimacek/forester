/*
 * Copyright (C) 2014 Martin Hruska
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

#include "vata_adapter.hh"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>

#include "connection_graph.hh"

VATAAdapter::VATAAdapter(TreeAut aut) : vataAut_(aut)
{}

VATAAdapter::VATAAdapter() : vataAut_()
{}

VATAAdapter::VATAAdapter(const VATAAdapter& adapter) : vataAut_(adapter.vataAut_)
{}

VATAAdapter::VATAAdapter(VATAAdapter&& adapter) : vataAut_(adapter.vataAut_)
{}

VATAAdapter::~VATAAdapter()
{}

VATAAdapter VATAAdapter::createTAWithSameTransitions(
    const VATAAdapter&         ta)
{
    //FA_DEBUG_AT(1,"Create TA with same transitions\n");
    return VATAAdapter(TreeAut(ta.vataAut_, false, false));
}

VATAAdapter* VATAAdapter::allocateTAWithSameTransitions(
    const VATAAdapter&         ta)
{
    //FA_DEBUG_AT(1,"Allocate TA with same transitions\n");
    return new VATAAdapter(TreeAut(ta.vataAut_, false, false));
}

VATAAdapter VATAAdapter::createTAWithSameFinalStates(
    const VATAAdapter&         ta,
    bool                       copyFinalStates)
{
    //FA_DEBUG_AT(1,"Create TA with same final states\n");
    return VATAAdapter(TreeAut(ta.vataAut_, true, copyFinalStates));
}

VATAAdapter* VATAAdapter::allocateTAWithSameFinalStates(
    const VATAAdapter&        ta,
    bool                      copyFinalStates)
{
    //FA_DEBUG_AT(1,"Allocate TA with same final states\n");
    return new VATAAdapter(TreeAut(ta.vataAut_, true, copyFinalStates));
}

VATAAdapter::iterator VATAAdapter::begin() const
{
    //FA_DEBUG_AT(1,"TA begin\n");
    return vataAut_.begin();
}

VATAAdapter::iterator VATAAdapter::end() const
{
    //FA_DEBUG_AT(1,"TA end\n");
	return vataAut_.end();
}

typename VATAAdapter::DownAccessor::Iterator VATAAdapter::begin(
        const size_t parent) const
{
    //FA_DEBUG_AT(1,"TA Down begin\n");
    return vataAut_[parent].begin();
}

typename VATAAdapter::DownAccessor::Iterator VATAAdapter::end(
        const size_t parent) const
{
    //FA_DEBUG_AT(1,"TA Down end\n");
    return vataAut_[parent].end();
}

typename VATAAdapter::AcceptTrans::Iterator VATAAdapter::accBegin() const
{
    //FA_DEBUG_AT(1,"TA Acc begin\n");
   return vataAut_.GetAcceptTrans().begin();
}

typename VATAAdapter::AcceptTrans::Iterator VATAAdapter::accEnd() const
{
    //FA_DEBUG_AT(1,"TA Acc end\n");
   return vataAut_.GetAcceptTrans().end();
}

VATAAdapter& VATAAdapter::operator=(const VATAAdapter& rhs)
{
    //FA_DEBUG_AT(1,"TA =\n");

		if (this != &rhs)
		{
			this->vataAut_ = rhs.vataAut_;
		}

    return *this;
}

void VATAAdapter::addTransition(
		const std::vector<size_t>&          children,
		const SymbolType&                   symbol,
		const size_t                        parent)
{
    //FA_DEBUG_AT(1,"TA add transition\n");

    this->vataAut_.AddTransition(children, symbol, parent);
}

void VATAAdapter::addTransition(const Transition& transition)
{
    //FA_DEBUG_AT(1,"TA add transition 1\n");
    this->vataAut_.AddTransition(transition);
}

const VATAAdapter::Transition VATAAdapter::getTransition(
    const std::vector<size_t>&          children,
		const SymbolType&                   symbol,
		const size_t                        parent)

{
    //FA_DEBUG_AT(1,"TA get transition\n");
    if (vataAut_.ContainsTransition(children, symbol, parent))
    {
        return Transition(parent, symbol, children);
    }

    assert(false);
    return Transition();
}

const label_type VATAAdapter::GetSymbol(const Transition& t)
{
    //FA_DEBUG_AT(1,"TA get symbol\n");
    return label_type(t.GetSymbol());
}

void VATAAdapter::addFinalState(const size_t state)
{
    //FA_DEBUG_AT(1,"TA add final state\n");
    vataAut_.SetStateFinal(state);
}

void VATAAdapter::addFinalStates(const std::set<size_t>& states)
{
    //FA_DEBUG_AT(1,"TA add final states\n");
    vataAut_.SetStatesFinal(states);
}

void VATAAdapter::addFinalStates(const std::unordered_set<size_t>& states)
{
    //FA_DEBUG_AT(1,"TA add final states 1\n");
    for (const size_t state : states)
    {
        vataAut_.SetStateFinal(state);
    }
}

bool VATAAdapter::isFinalState(const size_t state) const
{
    //FA_DEBUG_AT(1,"TA is final state\n");
    return vataAut_.IsStateFinal(state);
}

const std::unordered_set<size_t>& VATAAdapter::getFinalStates() const
{
    //FA_DEBUG_AT(1,"TA get final states\n");
    return vataAut_.GetFinalStates();
}

size_t VATAAdapter::getFinalState() const
{
		//FA_DEBUG_AT(1,"TA get final state\n");

	const std::unordered_set<size_t>& finalStates = this->getFinalStates();
	assert(1 == finalStates.size());
    return *finalStates.begin();
}

void VATAAdapter::eraseFinalStates()
{
	this->vataAut_.EraseFinalStates();
}

const VATAAdapter::Transition VATAAdapter::getAcceptingTransition() const
{
    assert(vataAut_.GetAcceptTrans().begin() != vataAut_.GetAcceptTrans().end());
    assert(++(vataAut_.GetAcceptTrans().begin()) == vataAut_.GetAcceptTrans().end());
    //FA_DEBUG_AT(1,"TA get accepting transitions\n");
    return *(vataAut_.GetAcceptTrans().begin());
}

VATAAdapter& VATAAdapter::unreachableFree(VATAAdapter& dst) const
{
    //FA_DEBUG_AT(1,"TA unreachable\n");
    dst.vataAut_ = std::move(vataAut_.RemoveUnreachableStates());
    return dst;
}

VATAAdapter& VATAAdapter::uselessAndUnreachableFree(VATAAdapter& dst) const
{
    //FA_DEBUG_AT(1,"TA useless\n");
    dst.vataAut_ = std::move(vataAut_.RemoveUselessStates());
    dst.vataAut_ = std::move(dst.vataAut_.RemoveUnreachableStates());
    return dst;
}

VATAAdapter& VATAAdapter::disjointUnion(
		VATAAdapter&                      dst,
		const VATAAdapter&                src,
		bool                              addFinalStates)
{
    //FA_DEBUG_AT(1,"TA disjoint\n");
    if (addFinalStates)
    {
        dst.vataAut_ = std::move(TreeAut::UnionDisjointStates(
            src.vataAut_, dst.vataAut_));
    }
    else
    { // if it is not needed copy the final states
      // it is sufficient to copy just the transitions
        for (const Transition& t : src.vataAut_)
        {
            dst.addTransition(t);
        }
    }

    return dst;
}

VATAAdapter& VATAAdapter::minimized(VATAAdapter& dst) const
{
    //FA_DEBUG_AT(1,"TA minimized\n");
    dst.vataAut_ = std::move(vataAut_.Reduce());
    return dst;
}

bool VATAAdapter::areTransitionsEmpty()
{
    //FA_DEBUG_AT(1,"TA are transitions empty\n");
    return vataAut_.AreTransitionsEmpty();
}

VATAAdapter& VATAAdapter::copyTransitions(VATAAdapter& dst) const
{
    //FA_DEBUG_AT(1,"TA copy transitions\n");
    CopyAllFunctor copyAllFunctor;
    dst.vataAut_.CopyTransitionsFrom(vataAut_, copyAllFunctor);

	return dst;
}

VATAAdapter& VATAAdapter::copyNotAcceptingTransitions(
        VATAAdapter&                       dst,
        const VATAAdapter&                 ta) const
{
    //FA_DEBUG_AT(1,"TA copy not accepting transitions\n");
    CopyNonAcceptingFunctor copyFunctor(ta);
    dst.vataAut_.CopyTransitionsFrom(vataAut_, copyFunctor);
		return dst;
}
	
void VATAAdapter::clear()
{
    //FA_DEBUG_AT(1,"TA clear\n");
    vataAut_.Clear();
}

bool VATAAdapter::subseteq(const VATAAdapter& a, const VATAAdapter& b)
{
    //FA_DEBUG_AT(1,"TA subseteq\n");
   	return TreeAut::CheckInclusion(a.vataAut_, b.vataAut_);
}

VATAAdapter& VATAAdapter::unfoldAtRoot(
	VATAAdapter&                   dst,
	const size_t                   newState,
	const bool                     registerFinalState) const
{
    //FA_DEBUG_AT(1,"TA unfoldAtRoot\n");
    if (registerFinalState)
    {
        dst.addFinalState(newState);
    }

    for (const Transition& trans : vataAut_)
    {
        dst.addTransition(trans);
        if (isFinalState(trans.GetParent()))
        {
            dst.addTransition(trans.GetChildren(), trans.GetSymbol(), newState);
        }
    }

    return dst;
}


VATAAdapter& VATAAdapter::unfoldAtRoot(
    VATAAdapter&                                  dst,
    const std::unordered_map<size_t, size_t>&     statesTranslator,
    bool                                          registerFinalState) const
{
    //FA_DEBUG_AT(1,"TA unfoldAtRoot1\n");
    this->copyTransitions(dst);
    for (auto state : this->getFinalStates())
    {
        std::unordered_map<size_t, size_t>::const_iterator j = statesTranslator.find(state);
        assert(j != statesTranslator.end());

        for (const Transition& trans : vataAut_[state])
        {
            dst.addTransition(trans.GetChildren(), trans.GetSymbol(), j->second);
        }

        if (registerFinalState)
        {
            dst.addFinalState(j->second);
        }
    }

    return dst;
}

void VATAAdapter::buildStateIndex(Index<size_t>& index) const
{
    //FA_DEBUG_AT(1,"TA buildStateIndex\n");

    for (auto state : vataAut_.GetUsedStates())
    {
        index.add(state);
    }
}


std::unordered_set<size_t> VATAAdapter::getUsedStates() const
{
    return vataAut_.GetUsedStates();
}

VATAAdapter::TreeAut::AcceptTrans VATAAdapter::getEmptyRootTransitions() const
{
    //FA_DEBUG_AT(1,"TA get empty root\n");
    assert(vataAut_.IsStateFinal(cEmptyRootTransIndex)
            && vataAut_.GetFinalStates().size() == 1);
    return vataAut_.GetAcceptTrans();
}

void VATAAdapter::copyReachableTransitionsFromRoot(
    const VATAAdapter&        src,
    const size_t&             rootState)
{
    //FA_DEBUG_AT(1,"TA copy reachable transitions from root\n");
    std::vector<size_t> stack;
    std::unordered_set<size_t> visited;

		stack.push_back(rootState);

    while(!stack.empty())
    {
        const size_t state = stack.back();
        stack.pop_back();
        visited.insert(state);
        for (const auto& t : src.vataAut_[state])
        {
            this->addTransition(t);
            for (const auto& ch : t.GetChildren())
            {
                if (visited.count(ch) == 0)
                {
                    stack.push_back(ch);
                }
            }
        }
    }
}

VATAAdapter& VATAAdapter::collapsed(
    VATAAdapter&                              dst,
    const std::unordered_map<size_t, size_t>& rel) const
{
    dst.vataAut_ = std::move(vataAut_.CollapseStates(rel));

    return dst;
}

VATAAdapter VATAAdapter::intersectionBU(
		const VATAAdapter&                   lhs,
		const VATAAdapter&                   rhs,
		VATA::AutBase::ProductTranslMap*     pTranslMap)
{
    return VATAAdapter(VATA::ExplicitTreeAut::IntersectionBU(lhs.vataAut_, rhs.vataAut_, pTranslMap));
}

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

std::ostream& operator<<(std::ostream& os, const VATAAdapter& ta)
{
    os << "TREE AUT " << std::endl;

    for (const auto t : ta.vataAut_)
    {
		if (ta.isFinalState(t.GetParent()))
		{
			os << "[" << stateToString(t.GetParent()) << "] " << VATAAdapter::GetSymbol(t) << " ";
		}
		else 
		{
			os << stateToString(t.GetParent()) << " " << VATAAdapter::GetSymbol(t) << " ";
		}
        for (auto s : t.GetChildren())
		{
			os << stateToString(s) << " ";
			//os << s << " " ;
		}
        os << "\n";
    }

    os.flush();
    return os;
}
