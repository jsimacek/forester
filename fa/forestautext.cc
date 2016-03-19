/*
 * Copyright (C) 2010 Jiri Simacek
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
#include "forestautext.hh"
#include "streams.hh"
#include "vata_adapter.hh"

bool FAE::subseteq(const FAE& lhs, const FAE& rhs)
{
	if (lhs.getRootCount() != rhs.getRootCount())
		return false;

	if (lhs.connectionGraph.data != rhs.connectionGraph.data)
		return false;

	for (size_t i = 0; i < lhs.getRootCount(); ++i)
	{
		if (!TreeAut::subseteq(*lhs.getRoot(i), *rhs.getRoot(i)))
			return false;
	}

	return true;
}


TreeAut& FAE::relabelReferences(
	TreeAut&                      dst,
	const TreeAut&                src,
	const std::vector<size_t>&    index)
{
	dst.addFinalStates(src.getFinalStates());
	for (const Transition& tr : src)
	{
		if (TreeAut::GetSymbol(tr)->isData())
			continue;

		std::vector<size_t> lhs;
		for (size_t state : tr.GetChildren())
		{
			const Data* data;
			if (this->isData(state, data))
			{
				if (data->isRef())
				{
					if (index[data->d_ref.root] != static_cast<size_t>(-1))
					{
						lhs.push_back(this->addData(dst, Data::createRef(index[data->d_ref.root], data->d_ref.displ)));
					}
					else
					{
						lhs.push_back(this->addData(dst, Data::createUndef()));
					}
				} else {
					lhs.push_back(this->addData(dst, *data));
				}
			} else
			{
				lhs.push_back(state);
			}
		}

		dst.addTransition(lhs, TreeAut::GetSymbol(tr), tr.GetParent());
	}

	return dst;
}


TreeAut& FAE::invalidateReference(
	TreeAut&                       dst,
	const TreeAut&                 src,
	size_t                         root)
{
	dst.addFinalStates(src.getFinalStates());
	for (const Transition& tr : src)
	{
		std::vector<size_t> lhs;
		for (size_t state : tr.GetChildren())
		{
			const Data* data;
			if (FAE::isData(state, data) && data->isRef(root))
			{
				lhs.push_back(this->addData(dst, Data::createUndef()));
			} else {
				lhs.push_back(state);
			}
		}
		if (!FAE::isRef(TreeAut::GetSymbol(tr), root))
			dst.addTransition(lhs, TreeAut::GetSymbol(tr), tr.GetParent());
	}
	return dst;
}


// TODO: I really need refactore 
// TODO: Error: You have to relable all states that
// were changed because of new data value
void FAE::setLabelsToValue(
	const size_t           root,
	const int              value,
	const size_t           bytesToSet)
{
	const auto& ta = this->getRoot(root);
	std::unordered_set<size_t> underRoot; 

	for (auto transIter = ta->accBegin(); transIter != ta->accEnd(); ++transIter)
	{
		const auto children = (*transIter).GetChildren();
		underRoot.insert(children.begin(), children.end());
	}

	class CopyFunctor : public VATA::ExplicitTreeAut::AbstractCopyF
	{
	private:
		const std::unordered_set<size_t>& underRoot_;
		const std::unordered_set<size_t>& finalStates_;

    public:
        CopyFunctor(
				const std::unordered_set<size_t>& underRoot,
				const std::unordered_set<size_t>& finalStates) :
			underRoot_ (underRoot), finalStates_(finalStates)
        {}

        bool operator()(const Transition& t)
        {
			const auto& label = *(TreeAut::GetSymbol(t));

			if ((underRoot_.count(t.GetParent()) && label.isData()) || finalStates_.count(t.GetParent()))
			{
				return false;
			}
			
			return true;
        }

	};

	// Copy transitions except the accepting ones and the ones under root
	std::shared_ptr<TreeAut> newTa = std::shared_ptr<TreeAut>(new TreeAut());
	ta->copyTransitionsWithFunctor(*newTa, CopyFunctor(underRoot, ta->getFinalStates()));
	newTa->addFinalStates(ta->getFinalStates());

	size_t setBytes = 0;
	for (auto transIter = ta->accBegin(); transIter != ta->accEnd(); ++transIter)
	{
		const auto& sels = *(TreeAut::GetSymbol(*transIter)->getSels());
		assert(sels.size() == (*transIter).GetChildren().size() || sels.size() - 1 == (*transIter).GetChildren().size());

		std::vector<size_t> children;
		size_t i = 0;
		for (const SelData& sel : sels)
		{
			auto childState = (*transIter).GetChildren().at(i);

			for (auto childIter = ta->begin(childState); childIter != ta->end(childState); ++childIter)
			{
				const auto& childLabel = *(TreeAut::GetSymbol(*childIter));
				if (!childLabel.isData())
				{
					continue; // wrong transition
				}
				else if (setBytes >= bytesToSet)
				{
					newTa->addTransition((*childIter).GetChildren(), TreeAut::GetSymbol(*childIter), childState);
				}
				else
				{
					childState = this->addData(*newTa, Data::createNewValue(childLabel.getData(), value));
				}
			}
			
			children.push_back(childState);
			++i;
			setBytes += sel.size;
		}
		newTa->addTransition(children, TreeAut::GetSymbol(*transIter), (*transIter).GetParent());
	}

	this->setRoot(root, newTa);

	return;
}


void FAE::freePosition(size_t root, const std::unordered_set<size_t>& rootsReferencedByVar)
{
	// Preconditions
	assert(root < roots_.size());
	assert(nullptr != roots_[root]);

	size_t pos;
	for (pos = 0; pos < roots_.size(); ++pos)
	{	// try to find a gap in the forest
		if (nullptr == roots_[pos] && !rootsReferencedByVar.count(pos))
		{
			break;
		}
	}

	if (roots_.size() == pos)
	{ // in case no gap was found
		roots_.push_back(std::shared_ptr<TreeAut>());
	}

	while (rootsReferencedByVar.count(pos))
	{
		roots_.push_back(std::shared_ptr<TreeAut>());
		++pos;
	}

	// 'pos' is now the target index
	assert(nullptr == roots_[pos]);

	std::vector<size_t> index(roots_.size());
	for (size_t i = 0; i < index.size(); ++i)
	{	// create the index
		index[i] = i;
	}

	index[root] = pos;
	index[pos] = static_cast<size_t>(-1);

	std::ostringstream os;
	utils::printCont(os, index);
	FA_DEBUG_AT(1,"Index: " << os.str());

	// relabel references to the TA that was at 'pos'
	for (size_t i = 0; i < roots_.size(); ++i)
	{
		if (nullptr != roots_[i])
		{	// if there is some TA
			roots_[i] = std::shared_ptr<TreeAut>(
				this->relabelReferences(roots_[i].get(), index));
		}
	}

	// swap the TA
	roots_[pos].swap(roots_[root]);

	FA_DEBUG_AT(1,"TA " << root << " moved to " << pos << ".");
	FA_DEBUG_AT(1,"FAE now: " << *this);
}

