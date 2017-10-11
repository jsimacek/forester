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

// DO NOT MOVE LINES IN THIS PART  --- BEGIN

// Forester headers
#include "config.h"
#include "streams.hh"

/**
 * @brief  Reports an error
 *
 * @param[in]  errMsg  The error message
 */
void reportErrorNoLocation(const char* errMsg)
{
	FA_ERROR(errMsg);
}

// DO NOT MOVE LINES IN THIS PART  --- END


// Standard library headers
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <set>
#include <algorithm>

// Code Listener headers
#include <cl/cl_msg.hh>
#include <cl/cldebug.hh>
#include <cl/clutil.hh>
#include <cl/code_listener.h>
#include <cl/storage.hh>
#include "../cl/ssd.h"

// Forester headers
#include "backward_run.hh"
#include "executionmanager.hh"
#include "fixpoint.hh"
#include "fixpointinstruction.hh"
#include "forestautext.hh"
#include "memplot.hh"
#include "microcode.hh"
#include "programconfig.hh"
#include "programerror.hh"
#include "restart_request.hh"
#include "symctx.hh"
#include "symexec.hh"
#include "svtrace_lite.hh"

using namespace ssd;

#if 0
void dumpOperandTypes(std::ostream& os, const cl_operand* op) {
	os << "operand:" << std::endl;
	cltToStream(os, op->type, false);
	os << "accessors:" << std::endl;
	const cl_accessor* acc = op->accessor;
	while (acc) {
		cltToStream(os, acc->type, false);
		acc = acc->next;
	}
}
#endif

// anonymous namespace
namespace
{
/**
 * @brief  Prints the trace to output stream
 *
 * @param[in,out]  os     The output stream
 * @param[in]      trace  The trace to be printed
 *
 * @returns  Modified stream
 */
std::ostream& printTrace(
	std::ostream&             os,
	const SymState::Trace&    trace)
{
	const CodeStorage::Insn* lastInsn = nullptr;

	for (auto it = trace.crbegin(); it != trace.crend(); ++it)
	{	// traverse in the reverse order
		const SymState& state = **it;

		assert(state.GetInstr());
		const AbstractInstruction& instr = *state.GetInstr();

		const CodeStorage::Insn* origInsn = instr.insn();
		if ((nullptr != origInsn) && (lastInsn != origInsn))
		{
			std::ostringstream oss;
			oss << *origInsn;

			std::string filename = MemPlotter::plotHeap(state, "trace", &origInsn->loc);
			lastInsn = origInsn;
			os << std::setw(50) << std::left
				<< Compiler::Assembly::insnToString(*origInsn) << " " << origInsn->loc.line << " // "
				<< origInsn->loc.file << ":" << std::setw(4) << std::left
				<< origInsn->loc.line << "|  " << filename << "\n";
		}
	}

	return os;
}


/**
 * @brief  Prints the microcode trace to output stream
 *
 * @param[in,out]  os     The output stream
 * @param[in]      trace  The trace to be printed
 *
 * @returns  Modified stream
 */
std::ostream& printUcodeTrace(
	std::ostream&             os,
	const SymState::Trace&    trace)
{
	const CodeStorage::Insn* lastInsn = nullptr;

	for (auto it = trace.crbegin(); it != trace.crend(); ++it)
	{	// traverse in the reverse order
		const SymState& state = **it;

		assert(state.GetInstr());
		const AbstractInstruction& instr = *state.GetInstr();
		const CodeStorage::Insn* clInsn = instr.insn();

		os << &state;

		os << std::setw(18);
		if (instr.isTarget())
		{
			std::ostringstream addrStream;
			addrStream << &instr;

			if ((nullptr != clInsn) && (clInsn != lastInsn)
				&& (clInsn->bb->front() == clInsn))
			{
				addrStream << " (" << clInsn->bb->name() << ")";
			}

			addrStream << ":";

			os << std::left << addrStream.str();
		}
		else
		{
			os << "";
		}

		std::ostringstream osInstr;
		osInstr << instr;

		os << std::setw(40) << std::left << osInstr.str();

		if ((nullptr != clInsn) && (lastInsn != clInsn))
		{
			lastInsn = clInsn;
			os << "; " << clInsn->loc.line << ": " << *clInsn;
		}

		os << "\n";
		//MemPlotter::plotHeap(state);
	}

	return os;
}

void dumpTrace(const SymState& state)
{
	std::vector<const SymState*> states;

	for (const SymState* pState = &state; pState != nullptr; pState = static_cast<const SymState*>(pState->GetParent()))
		states.push_back(pState);

	for (int i = states.size() - 1; i >= 0; --i)
	{
		FA_DEBUG_AT(2, "TRACE " << states.size() - i);

		if (nullptr != states[i]->GetInstr()->insn())
			FA_DEBUG_AT(2, states[i]->GetInstr()->insn()->loc << *states[i]->GetInstr()->insn());

		FA_DEBUG_AT(2, *states[i]);
	}
}

} // namespace


class SymExec::Engine
{
private:  // data members

	TreeAut ta_;
	TreeAut fixpoint_;
	BoxMan boxMan_;

	Compiler compiler_;
	Compiler::Assembly assembly_;

	ExecutionManager execMan_;

	const ProgramConfig &conf_;

	volatile bool dbgFlag_;
	volatile bool userRequestFlag_;

	std::vector<std::shared_ptr<const TreeAut>> predicates_;

protected:

	/**
	 * @brief  Prints boxes
	 *
	 * Method that prints all boxes from the box manager.
	 */
	void printBoxes() const
	{
		std::vector<const Box *> boxes;

		boxMan_.boxDatabase().asVector(boxes);

		std::map<std::string, const Box *> orderedBoxes;

		// reorder according to the name
		for (auto &box : boxes)
		{
			std::stringstream ss;

			ss << *static_cast<const AbstractBox *>(box);

			orderedBoxes.insert(std::make_pair(ss.str(), box));
		}

		for (auto &nameBoxPair : orderedBoxes)
		{
			FA_DEBUG_AT(1, nameBoxPair.first << ':' << std::endl << *nameBoxPair.second);
		}
	}

	/**
	 * @brief  Clears all fixpoints
	 */
	void clearFixpoints()
	{
		// clear all fixpoints
		for (auto instr : assembly_.code_)
		{
			if (instr->getType() == fi_type_e::fiFix)
			{
				// clear the fixpoint
				static_cast<FixpointInstruction *>(instr)->clear();
			}
		}
	}

	void clearFixpointsReverseStructures()
	{
		for (auto instr : assembly_.code_)
		{
			if (instr->getType() == fi_type_e::fiFix)
			{
				static_cast<FixpointInstruction *>(instr)->clearReverse();
			}
		}
	}


	void printInstructionInfo(const CodeStorage::Insn *insn, const SymState *state)
	{
		if (nullptr != insn)
		{    // in case current instruction IS an instruction
			FA_DEBUG_AT(2, SSD_INLINE_COLOR(C_LIGHT_RED, insn->loc << *insn));
			FA_DEBUG_AT(2, *state);
		}
		else
		{
			FA_DEBUG_AT(3, *state);
		}
	}

	void printRefinementInfo(const SymState *failPoint)
	{
		FA_DEBUG_AT(1, "The counterexample IS (PROBABLY) spurious");

		FA_DEBUG_AT(1, "Failing instruction: " << *failPoint->GetInstr());
		for (const auto &p : predicates_)
		{
			FA_DEBUG_AT(1, "Learnt predicate: " << *p);
		}

		FA_DEBUG_AT(1, "Predicates printed");
	}

	void printTraceInternal(
			const std::vector<const CodeStorage::Insn *> trace)
	{
		// prepare output
		std::streambuf *buf;
		std::ofstream of;
		if (conf_.traceFile.length() > 0)
		{
			of.open(conf_.traceFile, std::ofstream::out);
			buf = of.rdbuf();
		}
		else
		{
			buf = std::cerr.rdbuf();
		}

		std::ostream out(buf);
		SVTraceLite svPrinter;
		svPrinter.printTrace(trace, out);

		if (conf_.traceFile.length() > 0)
		{
			of.close();
		}
	}

	void printTrace(const cl_loc *location, const SymState::Trace &origTrace)
	{
		if (conf_.printSVTrace)
		{
			FA_LOG_MSG(location, "Printing SV-Comp trace");

			// preprocess trace
			std::vector<const CodeStorage::Insn *> trace;
			for (const auto &state : origTrace)
			{
				const CodeStorage::Insn *s = state->GetInstr()->insn();
				if (s != NULL)
				{
					trace.insert(trace.begin() + 0, s);
				}
			}

			printTraceInternal(trace);
		}

		if (conf_.printUcodeTrace)
		{
			FA_LOG_MSG(location, "Printing microcode trace");

			std::ostringstream oss;
			printUcodeTrace(oss, origTrace);
			Streams::traceUcode(oss.str().c_str());
		}
	}

	void reportRealError(const CodeStorage::Insn *insn, const ProgramError &e)
	{
		FA_NOTE("The counterexample IS real");
		if (nullptr != insn)
			FA_NOTE_MSG(&insn->loc, SSD_INLINE_COLOR(C_LIGHT_RED, *insn));
		if (nullptr != e.location())
			FA_ERROR_MSG(e.location(), e.what());
		else
			reportErrorNoLocation(e.what());
	}

	bool isAbstractionInstruction(AbstractInstruction *instr)
	{
		return fi_type_e::fiFix == instr->getType();
	}


	void addNewPredicates()
	{
		FA_DEBUG_AT(1, "Running the analysis with the folowing predicates:");
		for (AbstractInstruction *instr : this->GetAssembly().code_)
		{
			if (isAbstractionInstruction(instr))
			{
				FI_abs *absInstr = dynamic_cast<FI_abs *>(instr);
				if (nullptr != absInstr)
				{
					absInstr->addPredicate(predicates_);
					absInstr->printDebugInfoAboutPredicates();
				}
			}
		}
		FA_DEBUG_AT(1, "\n---------------------END---------------------------");
	}

	bool isSpurious(
			std::vector<std::shared_ptr<const TreeAut>> &predicate,
			const SymState::Trace &origTrace,
			SymState *&failPoint)
	{
		// check whether the counterexample is spurious and in case it is collect
		// some perhaps helpful information (failpoint and predicate)
		BackwardRun bwdRun(execMan_);

		return bwdRun.isSpuriousCE(origTrace, failPoint, predicate);
	}

	bool isSpurious(
			const SymState::Trace &origTrace)
	{
		std::vector<std::shared_ptr<const TreeAut>> predicate;
		SymState *failPoint = nullptr;

		return isSpurious(predicate, origTrace, failPoint);
	}

	void runBackwardRunAlone(const CodeStorage::Insn *insn, const ProgramError &e)
	{
		FA_LOG("Executing backward run...");

		if (isSpurious(e.state()->getTrace()))
		{
			FA_NOTE("The counterexample IS spurious");
		}
		else
		{
			reportRealError(insn, e);
		}
	}

	void assertAbstractionInstruction(AbstractInstruction *insn)
	{
		FI_abs *absInstr = dynamic_cast<FI_abs *>(insn);
		if (nullptr == absInstr)
		{
			assert(false);
		}
	}

	void symbolicExecutionRun()
	{
		SymState* state = nullptr;
		while (nullptr != (state = execMan_.dequeueDFS()))
        {	// process all states in the DFS order
            assert(nullptr != state);

            printInstructionInfo(state->GetInstr()->insn(), state);

            if (testAndClearUserRequestFlag())
            {
                FA_NOTE("Executed " << std::setw(7) << execMan_.statesEvaluated()
                    << " states and " << std::setw(7) << execMan_.pathsEvaluated()
                    << " paths so far.");
            }

            // run the state
            execMan_.execute(*state);
        }
	}

	bool processProgramError(const ProgramError& e)
	{
		bool shouldRefineAndContinue = false;

		const CodeStorage::Insn* insn = e.state()->GetInstr()->insn();
        if (nullptr != insn)
        {
            FA_DEBUG_AT(2, std::endl << *(e.state()->GetFAE()));
        }

        printTrace(e.location(), e.state()->getTrace());

        if (FA_BACKWARD_RUN && !FA_USE_PREDICATE_ABSTRACTION)
        {
            runBackwardRunAlone(insn, e);
        }
        else if (FA_USE_PREDICATE_ABSTRACTION)
        {	// in case we are using predicate abstraction
            FA_LOG("Executing backward run because of " << e.what());

            SymState* failPoint = nullptr;
            predicates_.clear();
            if (isSpurious(predicates_, e.state()->getTrace(), failPoint))
            {
                assert(!predicates_.empty());
                assert(nullptr != failPoint);
                assert(nullptr != failPoint->GetInstr());

                printRefinementInfo(failPoint);

                assertAbstractionInstruction(failPoint->GetInstr());
                addNewPredicates();

                clearFixpoints();
				clearFixpointsReverseStructures();

                shouldRefineAndContinue = true;
            }
        }

		if (!shouldRefineAndContinue)
		{
			reportRealError(insn, e);
		}

		dumpTrace(*e.state());

		return shouldRefineAndContinue;
	}


	/**
	 * @brief  The main execution loop
	 *
	 * This method is the main execution loop for the symbolic execution. It
	 * assumes that the microcode is already compiled, etc.
	 */
	bool mainLoop()
	{
		FA_DEBUG_AT(2, "creating empty heap ...");

		// create an empty heap
		std::shared_ptr<FAE> fae = std::shared_ptr<FAE>(
			new FAE(ta_, boxMan_));

		FA_DEBUG_AT(2, "scheduling initial state ...");

		// schedule the initial state for processing
		execMan_.init(
			DataArray(assembly_.regFileSize_, Data::createUndef()),
			fae,
			assembly_.code_.front()
		);

		try
		{	// expecting problems...
			symbolicExecutionRun();
			return true; // program is safe
		}
		catch (ProgramError& e)
		{
			assert(nullptr != e.state());

			const bool refineAndContinue = processProgramError(e);

			if (refineAndContinue)
			{
				return false; // run was unsuccessful
			}
			else
			{
				throw;
			}
		}
		catch (RestartRequest& e)
		{	// in case a restart is requested, clear all fixpoint computation points
			FA_NOTE("Restart requested - new box...");
			predicates_.clear();
			clearFixpoints();
			clearFixpointsReverseStructures();

			FA_DEBUG_AT(2, e.what());

			return false;
		}

		assert(false);
	}

public:   // methods

	/**
	 * @brief  The default constructor
	 *
	 * The default constructor.
	 */
	explicit Engine(const ProgramConfig& conf) :
		ta_{},
		fixpoint_{},
		boxMan_{},
		compiler_(fixpoint_, ta_, boxMan_),
		assembly_{},
		execMan_{},
		conf_(conf),
		dbgFlag_{false},
		userRequestFlag_{false},
		predicates_()
	{ }

	/**
	 * @brief  Loads types from a storage
	 *
	 * This method loads data types and function stackframes from the provided
	 * storage.
	 *
	 * @param[in]  stor  The code storage containing types
	 */
	void loadTypes(const CodeStorage::Storage& stor)
	{
		FA_DEBUG_AT(3, "loading types ...");

		// clear the box manager
		boxMan_.clear();

		// ************ infer data types' layouts ************
		for (const cl_type* type : stor.types)
		{	// for each data type in the storage
			std::vector<size_t> v;
			std::string name;
			if (type->name)
			{	// in case the structure has a name
				name = std::string(type->name);
			}
			else
			{	// in case the structure is nameless
				std::ostringstream ss;
				ss << type->uid;
				name = ss.str();
			}
			FA_DEBUG_AT(3, name);

			switch (type->code)
			{
				case cl_type_e::CL_TYPE_STRUCT: // for a structure
					NodeBuilder::buildNode(v, type);
					break;

				default: // for other types
					v.push_back(0);
					break;
			}

			boxMan_.createTypeInfo(name, v);
		}

		// ************ infer functions' stackframes ************
		for (const CodeStorage::Fnc* fnc : stor.fncs)
		{	// for each function in the storage, create a data structure representing
			// its stackframe
			std::vector<size_t> v;

			const SymCtx ctx(*fnc);
			for (const SelData& sel : ctx.GetStackFrameLayout())
			{	// create the stackframe
				v.push_back(sel.offset);
			}

			std::ostringstream ss;
			ss << "__@" << nameOf(*fnc) << ':' << uidOf(*fnc);

			FA_DEBUG_AT(3, ss.str());

			boxMan_.createTypeInfo(ss.str(), v);
		}

		// ************ compile layout of the block of global vars ************
		std::vector<const cl_type*> components;
		for (const CodeStorage::Var& var : stor.vars)
		{
			if (CodeStorage::EVar::VAR_GL == var.code)
			{
				components.push_back(var.type);
			}
		}

		std::vector<size_t> v;
		if (!components.empty())
		{ // in case the are some global variables
			NodeBuilder::buildNodes(v, components, 0);
		}
		else
		{	// in case there are no global variables, make one fake
			v.push_back(0);
		}

		boxMan_.createTypeInfo(GLOBAL_VARS_BLOCK_STR, v);
		FA_DEBUG_AT(1, "created box for global variables: "
			<< *boxMan_.getTypeInfo(GLOBAL_VARS_BLOCK_STR));
	}

#if 0
	void loadBoxes(const std::unordered_map<std::string, std::string>& db) {

		FA_DEBUG_AT(2, "loading boxes ...");

		for (auto p : db) {

			this->boxes.push_back((const Box*)boxMan_.loadBox(p.first, db));

			FA_DEBUG(p.first << ':' << std::endl << *(const FA*)this->boxes.back());

		}

		boxMan_.buildBoxHierarchy(this->hierarchy, this->basicBoxes);

	}
#endif

	void compile(const CodeStorage::Storage& stor, const CodeStorage::Fnc& entry)
	{
		compiler_.compile(assembly_, stor, entry);
	}

	const Compiler::Assembly& GetAssembly() const
	{
		return assembly_;
	}

	void run()
	{
		// Assertions
		assert(assembly_.code_.size());

		try
		{	// expect problems...
			while (!this->mainLoop())
			{	// while the analysis hasn't terminated
				FA_NOTE("Restarting the analysis...");
			}

			FA_NOTE("The program is SAFE.");

			// print out boxes
			this->printBoxes();

			for (auto instr : assembly_.code_)
			{
				if (conf_.printSVTrace && instr->getType() == fi_type_e::fiCheck)
				{
					if ((static_cast<FI_check *>(instr))->wasGarbageFound())
					{
						std::vector<const CodeStorage::Insn *> trace;
						printTraceInternal(trace);
					}
				}

				// print out all fixpoints
				if (instr->getType() != fi_type_e::fiFix)
				{
					continue;
				}

				if (instr->insn())
				{
					FA_DEBUG_AT(1, "fixpoint at " << instr->insn()->loc << std::endl
						<< (static_cast<FixpointInstruction*>(instr))->getFixPoint());
				} else
				{
					FA_DEBUG_AT(1, "fixpoint at unknown location" << std::endl
						<< (static_cast<FixpointInstruction*>(instr))->getFixPoint());
				}
			}

			// print out stats
			FA_DEBUG_AT(1, "forester has generated " << execMan_.statesEvaluated()
				<< " symbolic configuration(s) in " << execMan_.pathsEvaluated()
				<< " path(s) using " << boxMan_.boxDatabase().size() << " box(es)");
		}
		catch (const ProgramError& e)
		{ }
		catch (std::exception& e)
		{
			FA_DEBUG(e.what());

			this->printBoxes();

			throw;
		}
	}

	void run(const Compiler::Assembly& assembly)
	{
		assembly_ = assembly;

		try
		{
			this->run();
			assembly_.code_.clear();
		}
		catch (...)
		{
			assembly_.code_.clear();

			throw;
		}
	}

	void setDbgFlag()
	{
		dbgFlag_ = true;
	}

	void setUserRequestFlag()
	{
		userRequestFlag_ = true;
	}

	bool testAndClearUserRequestFlag()
	{
		bool oldValue = userRequestFlag_;
		userRequestFlag_ = false;
		return oldValue;
	}
};

SymExec::SymExec(const ProgramConfig& conf) :
	engine{new Engine(conf)}
{ }

SymExec::~SymExec()
{
	// Assertions
	assert(engine != nullptr);

	delete this->engine;
}

void SymExec::loadTypes(const CodeStorage::Storage& stor)
{
	// Assertions
	assert(engine != nullptr);

	this->engine->loadTypes(stor);
}

#if 0
void SymExec::loadBoxes(const std::unordered_map<std::string, std::string>& db) {
	this->engine->loadBoxes(db);
}
#endif

const Compiler::Assembly& SymExec::GetAssembly() const
{
	// Assertions
	assert(nullptr != engine);

	return this->engine->GetAssembly();
}

void SymExec::compile(const CodeStorage::Storage& stor,
	const CodeStorage::Fnc& main)
{
	// Assertions
	assert(engine != nullptr);

	this->engine->compile(stor, main);
}

void SymExec::run()
{
	// Assertions
	assert(engine != nullptr);

	this->engine->run();
}

void SymExec::run(const Compiler::Assembly& assembly)
{
	// Assertions
	assert(engine != nullptr);

	this->engine->run(assembly);
}

void SymExec::setDbgFlag()
{
	// Assertions
	assert(engine != nullptr);

	this->engine->setDbgFlag();
}

void SymExec::setUserRequestFlag()
{
	// Assertions
	assert(engine != nullptr);

	this->engine->setUserRequestFlag();
}
