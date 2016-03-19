#include "unfolding.hh"

#include "splitting.hh"

namespace
{
    TreeAut TARelabeledToBox(
        FAE&                           fae,
        const TreeAut&                 boxRoot,
        const std::vector<size_t>&     rootIndex)
    {
            TreeAut relabeledTA = fae.createTAWithSameBackend();
            fae.relabelReferences(relabeledTA, boxRoot, rootIndex);

            return relabeledTA;
    }

    TreeAut copyAndRelabelTAToBox(
        FAE&                           fae,
        const TreeAut&                 boxRoot,
        const std::vector<size_t>&     rootIndex)
    {
        TreeAut newTA = fae.createTAWithSameBackend();
//		this->fae.boxMan->adjustLeaves(tmp2, boxRoot);
        fae.unique(newTA, TARelabeledToBox(fae, boxRoot, rootIndex));

        return newTA;
    }

    template <class T>
    std::vector<T> vectorConcat(const std::vector<T>& v1, const std::vector<T>& v2)
    {
        auto res(v1);
        res.insert(res.end(), v2.begin(), v2.end());

        return res;
    }
}


void Unfolding::boxMerge(
		TreeAut&                       dst,
		const TreeAut&                 src,
		const TreeAut&                 boxRoot,
		const Box*                     box,
		const std::vector<size_t>&     rootIndex)
{
	TreeAut relabeledTA = copyAndRelabelTAToBox(this->fae, boxRoot, rootIndex);
    //		this->fae.boxMan->adjustLeaves(relabeledTA, boxRoot);
    // Copy the original and relabeled automaton to destination
    relabeledTA.copyNotAcceptingTransitions(dst, relabeledTA);
    dst.addFinalStates(relabeledTA.getFinalStates());
    src.copyNotAcceptingTransitions(dst, src);

    for (const size_t state : src.getFinalStates())
    {
        for (auto originalTrans = src.begin(state); originalTrans != src.end(state); ++originalTrans)
        {
            std::vector<size_t> boxLhs;
            std::vector<const AbstractBox*> boxLabel;
            getChildrenAndLabelFromBox(box, *originalTrans, boxLhs, boxLabel);

            for (auto relabeledTrans = relabeledTA.accBegin();
                 relabeledTrans != relabeledTA.accEnd(); ++relabeledTrans)
            {
                std::vector<size_t> newLhs =
                        vectorConcat(boxLhs, (*relabeledTrans).GetChildren());
                std::vector<const AbstractBox*> newLabel =
                        vectorConcat(boxLabel, TreeAut::GetSymbol(*relabeledTrans)->getNode());

                FA::reorderBoxes(newLabel, newLhs);

                if (nullptr == this->fae.boxMan->lookupLabel(newLabel).obj_)
                {
                    throw std::runtime_error("Cannot find box -- equivalent to invalid deref");
                }
                assert(nullptr != this->fae.boxMan->lookupLabel(newLabel).obj_);
                dst.addTransition(
                        newLhs,
                        static_cast<VATAAdapter::SymbolType>(
                                this->fae.boxMan->lookupLabel(newLabel)),
                        (*relabeledTrans).GetParent()
                );
            }
        }
    }
}

void Unfolding::getChildrenAndLabelFromBox(
    const Box*                         box,
    const TreeAut::Transition&         transition,
    std::vector<size_t>&               children,
    std::vector<const AbstractBox*>&   label)
{
    size_t childrenOffset = 0;
    if (box)
    {
        bool found = false;
        for (const AbstractBox* aBox : TreeAut::GetSymbol(transition)->getNode())
        {
            if (!aBox->isStructural())
            {
                label.push_back(aBox);
                continue;
            }
            const StructuralBox* b = static_cast<const StructuralBox*>(aBox);
            if (b != static_cast<const StructuralBox*>(box))
            {
                // this box is not interesting
                for (size_t k = 0; k < b->getArity(); ++k, ++childrenOffset)
                    children.push_back(transition.GetNthChildren(childrenOffset));
                label.push_back(b);
                continue;
            }
            childrenOffset += box->getArity();

            if (found)
                assert(false);

            found = true;
        }

        if (!found)
            assert(false);

    } else
    {
        children = transition.GetChildren();
        label = TreeAut::GetSymbol(transition)->getNode();
    }
}

void Unfolding::initRootRefIndex(
		std::vector<size_t>&          index,
		const Box*                    box,
		const TreeAut::Transition&    t)
{
    size_t lhsOffset = 0;

	// First initialize structures before unfolding
    for (const AbstractBox* aBox : TreeAut::GetSymbol(t)->getNode())
    {
        if (static_cast<const AbstractBox*>(box) != aBox)
        {
            lhsOffset += aBox->getArity();
        }
		else
		{
			for (size_t j = 0; j < box->getArity(); ++j)
			{
				const Data& data = this->fae.getData(t.GetNthChildren(lhsOffset + j));

				if (data.isUndef())
					index.push_back(static_cast<size_t>(-1));
				else
					index.push_back(data.d_ref.root);
			}

	        return;
		}
    }
}

void Unfolding::substituteOutputPorts(
		const std::vector<size_t>&    index,
		const size_t                  root,
		const Box*                    box)
{
    auto ta = std::shared_ptr<TreeAut>(this->fae.allocTA());

    this->boxMerge(*ta, *this->fae.getRoot(root), *box->getOutput(), box, index);

    this->fae.setRoot(root, ta);
    this->fae.connectionGraph.invalidate(root);
}

void Unfolding::substituteInputPorts(
		const std::vector<size_t>&    index,
		const Box*                    box)
{
    assert(box->getInputIndex() < index.size());
	const size_t aux = index.at(box->getInputIndex() + 1);

    assert(aux != static_cast<size_t>(-1));
    assert(aux < this->fae.getRootCount());

    TreeAut tmp = this->fae.createTAWithSameBackend();

    this->fae.getRoot(aux)->unfoldAtRoot(tmp, this->fae.freshState());
    this->fae.setRoot(aux, std::shared_ptr<TreeAut>(this->fae.allocTA()));

    this->boxMerge(*this->fae.getRoot(aux), tmp, *box->getInput(), nullptr, index);

    this->fae.connectionGraph.invalidate(aux);
}


bool Unfolding::updateBoxSignature(
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

	for (auto& boxOriginPair : boxSignature)
	{
		auto itBoolPair2 = prevSignature.insert(boxOriginPair);

		if (itBoolPair2.second)
		{
			changed = true;
			continue;
		}

		assert(
			boxOriginPair.second == itBoolPair2.first->second || boxOriginPair.second == static_cast<size_t>(-1) ||
			itBoolPair2.first->second == static_cast<size_t>(-1)
		);

		if (boxOriginPair.second == static_cast<size_t>(-1) && itBoolPair2.first->second != static_cast<size_t>(-1))
		{
			itBoolPair2.first->second = static_cast<size_t>(-1);
			changed = true;
		}
	}

	return changed;
}


void Unfolding::joinBoxSignature(BoxSignature& dst, const std::pair<const Box*, size_t>& boxOriginPair)
{
	auto itBoolPair = dst.insert(boxOriginPair);

	if (itBoolPair.second)
		return;

	// -1 indicates that box appears multiple times
	itBoolPair.first->second = static_cast<size_t>(-1);
}


void Unfolding::joinBoxSignature(BoxSignature& dst, const BoxSignature& src)
{
	for (auto& boxOriginPair : src)
		joinBoxSignature(dst, boxOriginPair);
}


bool Unfolding::processLhs(
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


void Unfolding::computeBoxSignatureMap(
		StateToBoxSignatureMap&	stateMap,
		const TreeAut&			ta)
{
	stateMap.clear();

	// the workset of transitions
	std::list<const TreeAut::Transition*> transitions;

	// compute the initial signatures for leaves, other signatures are cleared
	for (auto trans : ta)
	{	// traverse transitions of the TA
		if (TreeAut::GetSymbol(trans)->isData())
		{	// for data transitions
			stateMap.insert(std::make_pair(trans.GetParent(), BoxSignature()));
		} else
		{	// for non-data transitions
			transitions.push_back(&trans);
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

			if (!processLhs(boxSignature, t->GetChildren(), stateMap))
			{	// in case this transition cannot be processed because of some downward
				// states with missing box information
				continue;
			}

			for (const AbstractBox* box : TreeAut::GetSymbol(*t)->getNode())
			{	// for all boxes in the label
				assert(nullptr != box);

				if (box->isBox())
				{
					joinBoxSignature(boxSignature, std::make_pair(static_cast<const Box*>(box), t->GetParent()));
				}
			}

			changed |= updateBoxSignature(stateMap, t->GetParent(), boxSignature);
		}
	}
}

void Unfolding::unfoldBox(const size_t root, const Box* box)
{
    //std::cerr << "FA " << this->fae;
    //std::cerr << "Root Number " << root << std::endl;
    assert(root < this->fae.getRootCount());
    //std::cerr << "TA " << *(this->fae.getRoot(root));
    assert(nullptr != this->fae.getRoot(root));
    assert(nullptr != box);

    const TreeAut::Transition& t =
        this->fae.getRoot(root)->getAcceptingTransition();

    std::vector<size_t> index = { root };
	initRootRefIndex(index, box, t);

	substituteOutputPorts(index, root, box);

    if (!box->getInput())
    {
        return;
    }

	substituteInputPorts(index, box);
}

void Unfolding::unfoldBoxes(const size_t root, const std::set<const Box*>& boxes)
{
    for (std::set<const Box*>::const_iterator i = boxes.begin(); i != boxes.end(); ++i)
        this->unfoldBox(root, *i);
}

void Unfolding::unfoldSingletons(const size_t root)
{
    assert(root < this->fae.getRootCount());
    assert(nullptr != this->fae.getRoot(root));

	auto ta = *this->fae.getRoot(root);

	StateToBoxSignatureMap stateMap;

	computeBoxSignatureMap(stateMap, ta);

	std::unordered_map<size_t, std::set<const Box*>> stateToBoxSetMap;

	for (auto finalState : ta.getFinalStates())
	{
		auto boxSignatureIter = stateMap.find(finalState);

		assert(boxSignatureIter != stateMap.end());

		for (auto boxStatePair : boxSignatureIter->second)
		{
			if (boxStatePair.second != static_cast<size_t>(-1))
				stateToBoxSetMap.insert(std::make_pair(boxStatePair.second, std::set<const Box*>())).first->second.insert(boxStatePair.first);
		}
	}

	std::vector<size_t> splittingStates;
	for (auto stateBoxSetPair : stateToBoxSetMap)
		splittingStates.push_back(stateBoxSetPair.first);

	std::vector<size_t> roots;
	Splitting(this->fae).split(roots, root, splittingStates);

	auto i = 0;
	for (auto stateBoxSetPair : stateToBoxSetMap)
	{
		unfoldBoxes(roots[i], stateBoxSetPair.second);

		++i;
	}
}
