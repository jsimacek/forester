/*
 * Copyright (C) 2010 Kamil Dudka <kdudka@redhat.com>
 *
 * This file is part of predator.
 *
 * predator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * predator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with predator.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"
#include "symplot.hh"

#include <cl/cl_msg.hh>
#include <cl/storage.hh>
#include <cl/clutil.hh>

#include "symbt.hh"
#include "symheap.hh"
#include "symseg.hh"
#include "symutil.hh"
#include "util.hh"
#include "worklist.hh"

#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stack>
#include <string>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#ifndef DEBUG_SYMPLOT
#   define DEBUG_SYMPLOT 0
#endif

#ifndef SYMPLOT_STOP_AFTER_N_STATES
#   define SYMPLOT_STOP_AFTER_N_STATES 0
#endif

// singleton
class PlotEnumerator {
    public:
        static PlotEnumerator* instance() {
            return (inst_)
                ? (inst_)
                : (inst_ = new PlotEnumerator);
        }

        // generate kind of more unique name
        std::string decorate(std::string name);

    private:
        static PlotEnumerator *inst_;
        PlotEnumerator() { }
        // FIXME: should we care about the destruction?

    private:
        typedef std::map<std::string, int> TMap;
        TMap map_;
};

// /////////////////////////////////////////////////////////////////////////////
// implementation of PlotEnumerator
PlotEnumerator *PlotEnumerator::inst_ = 0;

std::string PlotEnumerator::decorate(std::string name) {
    // obtain a unique ID for the given name
    const int id = map_[name] ++;
#if SYMPLOT_STOP_AFTER_N_STATES
    if (SYMPLOT_STOP_AFTER_N_STATES < id) {
        CL_ERROR("SYMPLOT_STOP_AFTER_N_STATES (" << SYMPLOT_STOP_AFTER_N_STATES
                << ") exceeded, now stopping per user's request...");
        CL_TRAP;
    }
#endif

    // convert the ID to string
    std::ostringstream str;
    str << std::fixed
        << std::setfill('0')
        << std::setw(/* width of the ID suffix */ 4)
        << id;

    // merge name with ID
    name += "-";
    name += str.str();

#ifdef SYMPLOT_STOP_CONDITION
    if (SYMPLOT_STOP_CONDITION(name))
        CL_TRAP;
#endif

    return name;
}


// /////////////////////////////////////////////////////////////////////////////
// implementation of SymPlot
struct SymPlot::Private {
    const CodeStorage::Storage          *stor;
    const SymHeap                       *heap;
    std::ofstream                       dotStream;
    bool                                ok;
    const struct cl_loc                 *lw;
    WorkList<TValId>                  workList;
    std::set<TObjId>                    objDone;
    int                                 last;

    typedef std::pair<TObjId, TValId>                       TEdgeValueOf;
    std::vector<TEdgeValueOf>           evList;

    typedef std::pair<TValId, SymHeap::TOffVal>             TEdgeOffVal;
    std::vector<TEdgeOffVal>            ovList;

    typedef std::pair<TValId, TValId>                       TEdgeNeq;
    std::set<TEdgeNeq>                  neqSet;

    std::set<TObjId>                    heads;
    std::set<TObjId>                    nexts;
    std::set<TObjId>                    peers;

    bool openDotFile(const std::string &name);
    void closeDotFile();

    bool digFieldName(std::string &dst, TObjId obj);
    void plotNodeObj(TObjId obj);
    void plotNodeObjAnon(TObjId obj);
    void plotNodeValue(TValId val, enum cl_type_e code, const char *label);
    void plotNodeAux(int src, enum cl_type_e code, const char *label);
    void plotNeqZero(TValId val);

    void plotEdgePointsTo(TValId value, TObjId obj);
    void plotEdgeValueOf(TObjId obj, TValId value);
    void plotEdgeSub(TObjId obj, TObjId sub);

    void gobbleEdgeValueOf(TObjId obj, TValId value);
    void gobbleEdgeOffValue(TValId val, const SymHeap::TOffVal &ov);
    void gobbleEdgeNeq(TValId val1, TValId val2);
    void gobbleEdgeAlias(TValId val1, TValId val2);
    void emitPendingEdges();

    void plotSingleValue(TValId value);
    void plotZeroValue(TObjId obj);
    void digNext(TObjId obj);
    void openCluster(TObjId obj);

    bool handleCustomValue(TValId value);
    bool handleUnknownValue(TValId value);
    bool resolveValueOf(TValId *pDst, TObjId obj);
    bool resolvePointsTo(TObjId *pDst, TValId val);

    void digObjCore(TObjId obj);
    void digObj(TObjId obj);
    void digValues();
    void plotObj(TObjId obj);
    void plotCVar(CVar cVar);
};

#define SL_QUOTE(what) "\"" << what << "\""

bool SymPlot::Private::openDotFile(const std::string &plotName) {
    // compute a sort of unique file name
    PlotEnumerator *pe = PlotEnumerator::instance();
    std::string name(pe->decorate(plotName));
    std::string fileName(name + ".dot");

    // now please create the file
    this->dotStream.open(fileName.c_str(), std::ios::out);
    if (!this->dotStream) {
        CL_ERROR("unable to create file '" << fileName << "'");
        return false;
    }

    // open graph
    this->dotStream << "digraph " << SL_QUOTE(name) << " {" << std::endl
        << "\tlabel=<<FONT POINT-SIZE=\"18\">"
        << name << "</FONT>>;"                      << std::endl
        << "\tclusterrank=local;"                   << std::endl
        << "\tlabelloc=t;"                          << std::endl;

    CL_DEBUG("symplot: created dot file '" << fileName << "'");
    return this->dotStream;
}

void SymPlot::Private::closeDotFile() {
    // emit pending edges
    this->emitPendingEdges();

    // close graph
    this->dotStream << "}" << std::endl;
    if (!this->dotStream)
        this->ok = false;

    // close stream
    this->dotStream.close();
}

namespace {
    const char* prefixByCode(enum cl_type_e code) {
        switch (code) {
            case CL_TYPE_VOID:      return "void";
            case CL_TYPE_UNKNOWN:   return "?";
            case CL_TYPE_PTR:       return "*";
            case CL_TYPE_FNC:       return "fnc*";
            case CL_TYPE_STRUCT:    return "struct";
            case CL_TYPE_UNION:     return "union";
            case CL_TYPE_ARRAY:     return "array";
            case CL_TYPE_STRING:    return "string";
            case CL_TYPE_CHAR:      return "char";
            case CL_TYPE_BOOL:      return "bool";
            case CL_TYPE_INT:       return "int";
            case CL_TYPE_ENUM:      return "enum";
            default:                return "XXX";
        }
    }

    const char* colorByCode(enum cl_type_e code) {
        switch (code) {
            case CL_TYPE_VOID:      return "red";
            case CL_TYPE_UNKNOWN:   return "gray";
            case CL_TYPE_PTR:       return "blue";
            case CL_TYPE_FNC:       return "green";
            case CL_TYPE_STRUCT:    return "black";
            case CL_TYPE_UNION:     return "red";
            case CL_TYPE_ARRAY:     return "gray";
            case CL_TYPE_STRING:    return "gray";
            case CL_TYPE_CHAR:      return "gray";
            case CL_TYPE_BOOL:      return "gold";
            case CL_TYPE_INT:       return "gray";
            case CL_TYPE_ENUM:      return "gray";
            default:                return "black";
        }
    }
}

bool SymPlot::Private::digFieldName(std::string &dst, TObjId obj) {
    int nth;
    const TObjId parent = this->heap->objParent(obj, &nth);
    if (OBJ_INVALID == parent)
        // no chance since there is no parent
        return false;

    const struct cl_type *clt = this->heap->objType(parent);
    CL_BREAK_IF(!clt || !isComposite(clt));

    const char *name = clt->items[nth].name;
    if (!name)
        // anonymous unions involved?
        return false;

    dst = name;
    return true;
}

void SymPlot::Private::plotNodeObj(TObjId obj) {
    CL_BREAK_IF(obj <= 0);
    const struct cl_type *clt = this->heap->objType(obj);
    const enum cl_type_e code = (clt)
        ? clt->code
        : CL_TYPE_UNKNOWN;

    this->dotStream << "\t" << SL_QUOTE(obj);
    this->dotStream << " [shape=box";

    if (hasKey(this->heads, obj))
        this->dotStream << ", color=green, penwidth=3.0, style=dashed";
    else if (hasKey(this->peers, obj))
        this->dotStream << ", color=gold, penwidth=3.0, style=dashed";
    else if (hasKey(this->nexts, obj))
        this->dotStream << ", color=red, penwidth=3.0, style=dashed";
    else
        this->dotStream << ", color=" << colorByCode(code);

    // dig root object
    TObjId root = obj, next;
    while (OBJ_INVALID != (next = this->heap->objParent(root)))
        root = next;

    if (this->heap->cVar(0, root))
        // colorize on-stack object
        this->dotStream << ", fontcolor=blue";

    else
        // colorize heap (sub)object
        this->dotStream << ", fontcolor=red";

    this->dotStream << ", label=\"[" << prefixByCode(code) << "] #" << obj;

    CVar cVar;
    if (this->heap->cVar(&cVar, obj)) {
        this->dotStream << " (#" << cVar.uid;
        const CodeStorage::Var &var = this->stor->vars[cVar.uid];
        std::string name = var.name;
        if (!name.empty())
            this->dotStream << " - " << name;
        if (1 < cVar.inst)
            this->dotStream << ", inst = " << cVar.inst;
        this->dotStream << ")";
    }

    std::string filedName;
    if (digFieldName(filedName, obj))
        this->dotStream << " ." << filedName;

    this->dotStream << "\"];" << std::endl;
}

void SymPlot::Private::plotNodeObjAnon(TObjId obj) {
    const int size = this->heap->objSizeOfAnon(obj);
    CL_BREAK_IF(size <= 0);

    this->dotStream << "\t" << SL_QUOTE(obj)
        << " [shape=box, color=red, fontcolor=red"
        << ", label=\"RAW " << size << "B\"];"
        << std::endl;
}

void SymPlot::Private::plotNodeValue(TValId val, enum cl_type_e code,
                                     const char *label)
{
    // visualize the count of references as pen width
    const float pw = static_cast<float>(this->heap->usedByCount(val));

    this->dotStream << "\t" << SL_QUOTE(val)
        << " [shape=ellipse"
        << ", penwidth=" << pw
        << ", color=" << colorByCode(code)
        << ", fontcolor=green"
        << ", label=\"[" << prefixByCode(code) << "] #" << val;

    if (label)
        this->dotStream << " [" << label << "]";

    this->dotStream << "\"];" << std::endl;
}

void SymPlot::Private::plotNodeAux(int src, enum cl_type_e code,
                                   const char *label)
{
    const int id = ++(this->last);
    this->dotStream << "\t"
        << SL_QUOTE("lonely" << id)
        << " [shape=plaintext"
        << ", fontcolor=" << colorByCode(code)
        << ", label=" << SL_QUOTE(label) << "];"
        << std::endl;

    this->dotStream << "\t"
        << SL_QUOTE(src) << " -> " << SL_QUOTE("lonely" << id)
        << " [color=" << colorByCode(code)
        << "];" << std::endl;
}

void SymPlot::Private::plotNeqZero(TValId val) {
    const char *label = "!= 0";

    const int id = ++(this->last);
    this->dotStream << "\t"
        << SL_QUOTE("lonely" << id)
        << " [shape=plaintext"
        << ", fontcolor=" << colorByCode(CL_TYPE_VOID)
        << ", label=" << SL_QUOTE(label) << "];"
        << std::endl;

    this->dotStream << "\t"
        << SL_QUOTE(val) << " -> " << SL_QUOTE("lonely" << id)
        << " [color=" << colorByCode(CL_TYPE_BOOL)
        << "];" << std::endl;
}

void SymPlot::Private::plotEdgePointsTo(TValId value, TObjId obj) {
    this->dotStream << "\t" << SL_QUOTE(value) << " -> " << SL_QUOTE(obj)
        << " [color=green, fontcolor=green, label=\"pointsTo\"];"
        << std::endl;
}

void SymPlot::Private::plotEdgeValueOf(TObjId obj, TValId value) {
    this->dotStream << "\t" << SL_QUOTE(obj) << " -> " << SL_QUOTE(value)
        << " [color=blue, fontcolor=blue];"
        << std::endl;
}

void SymPlot::Private::plotEdgeSub(TObjId obj, TObjId sub) {
    this->dotStream << "\t" << SL_QUOTE(obj) << " -> " << SL_QUOTE(sub)
        << " [color=gray, style=dotted, arrowhead=open"
        << ", fontcolor=gray, label=\"field\"];"
        << std::endl;
}

void SymPlot::Private::gobbleEdgeValueOf(TObjId obj, TValId value) {
    TEdgeValueOf edge(obj, value);
    this->evList.push_back(edge);
}

void SymPlot::Private::gobbleEdgeOffValue(TValId val,
                                          const SymHeap::TOffVal &ov)
{
    TEdgeOffVal edge(val, ov);
    this->ovList.push_back(edge);
}

void SymPlot::Private::gobbleEdgeNeq(TValId val1, TValId val2) {
    // Neq predicates induce a symmetric relation, let's handle them such
    sortValues(val1, val2);

    TEdgeNeq edge(val1, val2);
    this->neqSet.insert(edge);
}

void SymPlot::Private::emitPendingEdges() {
    // plot all valueOf edges
    BOOST_FOREACH(const TEdgeValueOf &edge, this->evList) {
        this->plotEdgeValueOf(edge.first, edge.second);
    }

    // plot all off-value edges
    BOOST_FOREACH(const TEdgeOffVal &edge, this->ovList) {
        const TValId dst = edge.first;
        const SymHeap::TOffVal &ov = edge.second;
        this->dotStream << "\t" << SL_QUOTE(ov.first) << " -> " << SL_QUOTE(dst)
            << " [color=red, fontcolor=red,"
            << " label=\"[+" << ov.second << "]\"];"
            << std::endl;
    }

    // plot Neq edges
    std::set<TEdgeNeq> neqDone;
    BOOST_FOREACH(const TEdgeNeq &edge, this->neqSet) {
        const TValId v1 = edge.first;
        const TValId v2 = edge.second;
        if (!insertOnce(neqDone, TEdgeNeq(v1, v2)))
            continue;

        this->dotStream << "\t" << SL_QUOTE(v1) << " -> " << SL_QUOTE(v2)
            << " [color=gold, fontcolor=red, label=\"Neq\", arrowhead=none];"
            << std::endl;
    }

    // cleanup for next wheel
    this->evList.clear();
    this->ovList.clear();
    this->neqSet.clear();
}

void SymPlot::Private::plotSingleValue(TValId value) {
    if (value <= 0) {
        this->plotNodeValue(value, CL_TYPE_UNKNOWN, 0);
        return;
    }

    // visualize off-value relations
    SymHeap::TOffValCont offValues;
    this->heap->gatherOffValues(offValues, value);
    BOOST_FOREACH(const SymHeap::TOffVal &ov, offValues) {
        this->workList.schedule(ov.first);
        if (ov.second < 0)
            // we came to the predicate from the less interesting side; let's
            // just wait for the value on the opposite side of the predicate
            continue;

        this->gobbleEdgeOffValue(value, ov);
    }

    // traverse all Neq predicates
    TValList relatedVals;
    this->heap->gatherRelatedValues(relatedVals, value);
    BOOST_FOREACH(TValId peer, relatedVals) {
        if (0 < peer)
            this->workList.schedule(peer);

        if (this->heap->SymHeapCore::proveNeq(value, peer)) {
            if (VAL_NULL == peer) {
                // 'value' is said to be non-zero
                this->plotNeqZero(value);
                continue;
            }
            else if (VAL_NULL < peer) {
                // regular Neq predicate
                this->gobbleEdgeNeq(value, peer);
                continue;
            }
        }

        CL_WARN("SymPlot: unhandled predicate over values #"
                << value << " and #" << peer);
    }

    const struct cl_type *clt = 0;
    const TObjId target = this->heap->pointsTo(value);
    if (0 < target)
        clt = this->heap->objType(target);

    const enum cl_type_e code = (clt)
        ? clt->code
        : CL_TYPE_UNKNOWN;

    this->plotNodeValue(value, code, 0);
}

void SymPlot::Private::plotZeroValue(TObjId obj) {
    const struct cl_type *clt = this->heap->objType(obj);
    CL_BREAK_IF(!clt);

    const enum cl_type_e code = clt->code;
    switch (code) {
        case CL_TYPE_INT:
            this->plotNodeAux(obj, code, "[int] 0");
            break;

        case CL_TYPE_PTR:
            this->plotNodeAux(obj, code, "NULL");
            break;

        case CL_TYPE_BOOL:
            this->plotNodeAux(obj, code, "FALSE");
            break;

        default:
            this->plotNodeAux(obj, CL_TYPE_UNKNOWN, "[unknown type] 0");
    }
}

void SymPlot::Private::digNext(TObjId obj) {
    EObjKind kind = this->heap->objKind(obj);
    switch (kind) {
        case OK_CONCRETE:
        case OK_HEAD:
        case OK_PART:
            return;

        case OK_MAY_EXIST:
        case OK_SLS:
        case OK_DLS:
            break;
    }

    const BindingOff &off = this->heap->objBinding(obj);
    const TOffset offHead = off.head;
    if (offHead) {
        const TObjId objHead = compObjByOffset(*this->heap, obj, offHead);
        CL_BREAK_IF(objHead <= 0);

        // store 'head' pointer object
        this->heads.insert(objHead);
    }

    const TObjId objNext = ptrObjByOffset(*this->heap, obj, off.next);
    CL_BREAK_IF(objNext <= 0);

    // store 'next' poitner object
    this->nexts.insert(objNext);
    if (OK_DLS != kind)
        return;

    const TObjId objPeer = ptrObjByOffset(*this->heap, obj, off.prev);
    CL_BREAK_IF(objPeer <= 0);

    // store 'peer' pointer object
    this->peers.insert(objPeer);
}

void SymPlot::Private::openCluster(TObjId obj) {
    std::string label;
    if (this->heap->objIsProto(obj))
        label = "[prototype] ";

#ifndef NDEBUG
    const char *color, *pw;
#else
    const char *color = "", *pw = "";
#endif

    const struct cl_type *clt = this->heap->objType(obj);
    CL_BREAK_IF(!clt);

    const EObjKind kind = this->heap->objKind(obj);
    switch (kind) {
        case OK_CONCRETE:
        case OK_PART:
            color = (CL_TYPE_UNION == clt->code)
                ? "red"
                : "black";
            pw = "1.0";
            break;

        case OK_MAY_EXIST:
            label += "MAY_EXIST";
            color = "blue";
            pw = "3.0";
            break;

        case OK_HEAD:
            label += "head";
            color = "green";
            pw = "2.0";
            break;

        case OK_SLS:
            label += "SLS";
            color = "red";
            pw = "3.0";
            break;

        case OK_DLS:
            label += "DLS/2";
            color = "gold";
            pw = "3.0";
            break;
    }

    this->dotStream
        << "subgraph \"cluster" <<(++last)<< "\" {" << std::endl
        << "\trank=same;"                           << std::endl
        << "\tlabel=" << SL_QUOTE(label) << ";"     << std::endl
        << "\tcolor=" << color << ";"               << std::endl
        << "\tfontcolor=" << color << ";"           << std::endl
        << "\tbgcolor=gray98;"                      << std::endl;

    this->dotStream
        << "\tstyle=dashed;"                        << std::endl
        << "\tpenwidth=" << pw << ";"               << std::endl;
}

bool SymPlot::Private::handleCustomValue(TValId value) {
    using namespace CodeStorage;

    const int cVal = this->heap->valGetCustom(value);
    if (-1 == cVal)
        return false;

    const CodeStorage::FncDb &fncs = this->stor->fncs;
    const Fnc *fnc = fncs[cVal];
    CL_BREAK_IF(!fnc);

    const char *fncName = nameOf(*fnc);
    CL_BREAK_IF(!fncName);

    std::string name(fncName);
    name += "()";

    this->plotNodeValue(value, CL_TYPE_FNC, 0);
    this->plotNodeAux(value, CL_TYPE_FNC, name.c_str());
    return true;
}

bool SymPlot::Private::handleUnknownValue(TValId value) {
    const EUnknownValue code = this->heap->valGetUnknown(value);
    switch (code) {
        case UV_KNOWN:
        case UV_ABSTRACT:
            return false;

        case UV_UNINITIALIZED:
            this->plotNodeAux(value, CL_TYPE_UNKNOWN, "UV_UNINITIALIZED");
            return true;

        case UV_UNKNOWN:
            this->plotNodeAux(value, CL_TYPE_UNKNOWN, "UV_UNKNOWN");
            return true;

        case UV_DONT_CARE:
            this->plotNodeAux(value, CL_TYPE_UNKNOWN, "UV_DONT_CARE");
            return true;
    }

    CL_TRAP;
    return true;
}

bool SymPlot::Private::resolveValueOf(TValId *pDst, TObjId obj) {
    CL_BREAK_IF(obj < 0);

    // avoid duplicates
    if (hasKey(this->objDone, obj))
        return false;
    this->objDone.insert(obj);

    const TValId value = this->heap->valueOf(obj);
    switch (value) {
        case VAL_INVALID:
            this->plotNodeAux(obj, CL_TYPE_VOID, "VAL_INVALID");
            return false;

        case VAL_NULL /* = VAL_FALSE*/:
            this->plotZeroValue(obj);
            return false;

        case VAL_TRUE:
            this->plotNodeAux(obj, CL_TYPE_BOOL, "TRUE");
            return false;

        case VAL_DEREF_FAILED:
            this->plotNodeAux(obj, CL_TYPE_VOID, "UV_DEREF_FAILED");
            return false;

        default:
            break;
    }

    *pDst = value;
    return true;
}

bool SymPlot::Private::resolvePointsTo(TObjId *pDst, TValId value) {
    if (this->handleUnknownValue(value))
        return false;

    if (this->handleCustomValue(value))
        return false;

    const TObjId obj = this->heap->pointsTo(value);
    switch (obj) {
        case OBJ_INVALID:
            this->plotNodeAux(value, CL_TYPE_VOID, "INVALID");
            return false;

        case OBJ_DEREF_FAILED:
            this->plotNodeAux(value, CL_TYPE_VOID, "DEREF_FAILED");
            return false;

        case OBJ_DELETED:
            this->plotNodeAux(value, CL_TYPE_VOID, "DELETED");
            return false;

        case OBJ_LOST:
            this->plotNodeAux(value, CL_TYPE_VOID, "LOST");
            return false;

        case OBJ_UNKNOWN:
            this->plotNodeAux(value, CL_TYPE_UNKNOWN, "?");
            return false;

        case OBJ_RETURN:
        default:
            *pDst = obj;
            return true;
    }
}

class ObjectDigger {
    private:
        SymPlot::Private    *const self_;
        const TObjId        root_;
        unsigned            level_;

    public:
        ObjectDigger(SymPlot::Private *self, TObjId root):
            self_(self),
            root_(root),
            level_(0)
        {
            const struct cl_type *clt = self_->heap->objType(root);
            this->operate(TFieldIdxChain(), clt);
        }

        ~ObjectDigger() {
            // finally close all pending clusters
            this->setupNestLevel(0);
        }

        bool operator()(TFieldIdxChain ic, const struct cl_type_item *item) {
            this->operate(ic, item->type);
            return /* continue */ true;
        }

    private:
        void setupNestLevel(unsigned targetLevel) {
            for(; targetLevel < level_; --level_)
                self_->dotStream << "}" << std::endl;

            level_ = targetLevel;
        }

        void operate(TFieldIdxChain ic, const struct cl_type *clt);
};

void ObjectDigger::operate(TFieldIdxChain ic, const struct cl_type *clt) {
    CL_BREAK_IF(!clt);

    const SymHeap &sh = *self_->heap;
    const TObjId obj = subObjByChain(sh, root_, ic);
    CL_BREAK_IF(obj <= 0);

    // first close all pending clusters
    this->setupNestLevel(ic.size());

    if (!isComposite(clt)) {
        self_->plotNodeObj(obj);

        TValId value;
        if (self_->resolveValueOf(&value, obj)) {
            self_->gobbleEdgeValueOf(obj, value);
            self_->workList.schedule(value);
        }

        return;
    }

    // avoid duplicates
    self_->objDone.insert(obj);

    self_->digNext(obj);
    self_->openCluster(obj);
    self_->plotNodeObj(obj);
    for (int i = 0; i < clt->item_cnt; ++i) {
        const TObjId sub = self_->heap->subObj(obj, i);
        const struct cl_type *subClt = self_->heap->objType(sub);
        if (!subClt || (isComposite(subClt) && !subClt->item_cnt))
            // skip empty structures/unions
            continue;

        if (!hasKey(self_->objDone, sub))
            self_->plotEdgeSub(obj, sub);
    }
}

void SymPlot::Private::digObjCore(TObjId obj) {
    const struct cl_type *clt = this->heap->objType(obj);
    if (!clt) {
        this->plotNodeObjAnon(obj);
        return;
    }

    ObjectDigger visitor(this, obj);
    traverseTypeIc<TFieldIdxChain>(clt, visitor, /* digOnlyStructs */ true);
}

void SymPlot::Private::digObj(TObjId obj) {
    // seek root, in order to draw the whole object, even if the root is not
    // pointed from anywhere
    obj = objRoot(*this->heap, obj);

    if (OK_DLS != this->heap->objKind(obj)) {
        this->digObjCore(obj);
        return;
    }

    const TObjId peer = dlSegPeer(*this->heap, obj);
    const char *label = (this->heap->objIsProto(obj))
        ? "[prototype] DLS"
        : "DLS";

    // open a cluster
    this->dotStream
        << "subgraph \"cluster" <<(++last)<< "\" {" << std::endl
        << "\tlabel=" << SL_QUOTE(label)            << std::endl
        << "\tcolor=gold;"                          << std::endl
        << "\tfontcolor=gold;"                      << std::endl
        << "\tstyle=dashed;"                        << std::endl;

    // plot the two parts of a DLS into the cluster
    this->digObjCore(obj);
    if (0 < peer)
        this->digObjCore(peer);

    // close the cluster
    this->dotStream << "}" << std::endl;
}

void SymPlot::Private::digValues() {
    TValId value;
    while (workList.next(value)) {
        // plot the value itself
        this->plotSingleValue(value);

        if (value <= 0)
            // bare value can't be followed
            continue;

        TObjId obj = this->heap->valGetCompositeObj(value);
        if (OBJ_INVALID != obj) {
            // dig composite object and eventually schedule the values inside
            this->digObj(obj);
            continue;
        }

        // check the value inside
        if (!this->resolvePointsTo(&obj, value))
            // bare value can't be followed
            continue;

        // plot the pointing object and the corresponding "valueOf" edge
        this->plotNodeObj(obj);
        this->plotEdgePointsTo(value, obj);

        // follow values inside the object
        this->digObj(obj);
    }
}

void SymPlot::Private::plotObj(TObjId obj) {
    // plot the variable itself
    this->plotNodeObj(obj);

    // look for the value inside
    TValId value;
    if (!this->resolveValueOf(&value, obj))
        // we got a bare value, which can't be followed, so we're done
        return;

    if (OBJ_INVALID == this->heap->valGetCompositeObj(value)) {
        // connect the variable node with its value
        this->plotEdgeValueOf(obj, value);

        // dig the target value recursively and plot (if not already)
        this->workList.schedule(value);
    }
    else
        // dig composite object and eventually schedule the values inside
        this->digObj(obj);

    this->digValues();
}

void SymPlot::Private::plotCVar(CVar cVar) {
    // CodeStorage variable lookup
    const CodeStorage::Var &var = this->stor->vars[cVar.uid];
    this->lw = &var.loc;
#if DEBUG_SYMPLOT
    CL_DEBUG_MSG(this->lw, "-X- plotting stack variable: #" << var.uid
            << " (" << var.name << ")" );
#endif

    // SymbolicHeap variable lookup
    const TObjId obj = this->heap->objByCVar(cVar);
    if (OBJ_INVALID == obj)
        CL_DEBUG_MSG(this->lw, "objByCVar lookup failed");

    // plot as regular heap object
    this->plotObj(obj);
}

SymPlot::SymPlot(const CodeStorage::Storage &stor, const SymHeap &heap):
    d(new Private)
{
    d->stor = &stor;
    d->heap = &heap;
    d->last = 0;
}

SymPlot::~SymPlot() {
    delete d;
}

bool SymPlot::plot(const std::string &name) {
    // create dot file
    d->ok = true;
    if (!d->openDotFile(name))
        return false;

    // go through all program variables
    SymHeap::TContCVar cVars;
    d->heap->gatherCVars(cVars);
    BOOST_FOREACH(CVar cv, cVars) {
        d->plotCVar(cv);
    }

    // plot also all dangling objects, although we are not happy to see them
    TObjList roots;
    d->heap->gatherRootObjs(roots);
    BOOST_FOREACH(const TObjId obj, roots) {
        if (!hasKey(d->objDone, obj))
            d->plotObj(obj);
    }

    // close dot file
    d->closeDotFile();
    return d->ok;
}

bool SymPlot::plotHeapValue(const std::string &name, TValId value) {
    // create dot file
    d->ok = true;
    if (!d->openDotFile(name))
        return false;

    // plot by value
    d->workList.schedule(value);
    d->digValues();

    // close dot file
    d->closeDotFile();
    return d->ok;
}

bool SymPlot::plotStackFrame(const std::string           &name,
                             const CodeStorage::Fnc      &fnc,
                             const SymBackTrace          *bt)
{
    using namespace CodeStorage;

    // create dot file
    d->ok = true;
    if (!d->openDotFile(name))
        return false;

    d->lw = &fnc.def.loc;
#if DEBUG_SYMPLOT
    CL_DEBUG_MSG(d->lw, "-X- plotting stack frame of " << nameOf(fnc) << "():");
#endif

    // go through all stack variables
    BOOST_FOREACH(const int uid, fnc.vars) {
        const int nestLevel = bt->countOccurrencesOfFnc(uidOf(fnc));
        const CVar cVar(uid, nestLevel);
        d->plotCVar(cVar);
    }

    // close dot file
    d->closeDotFile();
    return d->ok;
}

// vim: tw=80
