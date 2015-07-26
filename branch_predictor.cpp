#include "bp.h"
#include <map>
#include<stdio.h>
#include <vector>
#include <math.h>

#define BUBBLE 0x5033


// Modify this line to substitute YOUR branch predictor for the
// default BTB...
#define BRANCH_PREDICTOR tlaBP
#define HIST_LEN 2
#define CONSTRAIN 170


//
// Custom Branch Predictors
//

//  TODO: add your BranchPredictor implementations here!

class PredictorTemplate : public BranchPredictor
{
private:
	// Define your predictor's internal data structures here
public:
	PredictorTemplate ( struct bp_io& io ) : BranchPredictor ( io )
	{
		// Initialize your predictor state here
	}

	~PredictorTemplate ( ) 
	{
		// delete any structures you heap-allocated here
	}

	uint32_t predict_fetch ( uint32_t pc )
	{
		// Predict the next PC here (or return 0 to predict not taken)
		return 0;
	}


	void update_execute ( uint32_t pc, 
	uint32_t pc_next, 
	bool mispredict,
	bool is_brjmp,
	uint32_t inst )
	{
		if ( is_brjmp )
		{
			// update your predictor state here
		}

		// Note that this function is called for every inst, not just br/jmp!
	}
};

  /******************************************/
 /****** Two-level Adaptive Predictor ******/
/******************************************/
class tlaBP : public BranchPredictor
{
private:

	struct Branch { // saturating counter and pattern frequency table
		uint32_t pc;
		uint32_t counter; // saturating counter
		std::map<uint32_t, bool> history; // taken history
		std::map<uint32_t, uint32_t> pattern; // pattern frequency table

		Branch(uint32_t _pc) {
			pc = _pc; // set pc
			counter = 0; // start counter at 0
			// The following creates a table entry for every possible
			// pattern that is HIST_LEN bits long.
			uint32_t maxPat = 0;
			for (uint32_t i = 0; i < HIST_LEN; i++) {
				maxPat += pow(2, i);
			}
			for (uint32_t i = 0; i <= maxPat; i++) {
				pattern[i] = 0; // set initial pattern frequencies to 0
			}
		}
	};
	
	typedef std::vector<Branch> branches;
	branches branchList; // stores all Branch structs
	std::map<uint32_t, uint32_t> table;
	bool taken;

public:
	tlaBP ( struct bp_io& _io ) : BranchPredictor ( _io ) { }
	~tlaBP ( ) {
		taken = 0;
	}

	uint32_t predict_fetch ( uint32_t pc ) {
		if (!exists(pc)) {
			return 0;
		}
		uint32_t prediction = predict(pc); // pattern prediction
		if (prediction == 2 || prediction == 3) { // predict taken
			return table[pc];
		}
		if (prediction == 0 || prediction == 1) { // predict not taken
			return 0;
		}
		return 0;
	}

	// Returns the two-level adaptive predictor counter.
	uint32_t predict(uint32_t pc) {
		uint32_t branchIndex = findBranch(pc);
		uint32_t patID = _toInt(branchList[branchIndex].history);
		return branchList[branchIndex].pattern[patID];
	}

	// Converts the boolean branch history into a unique integer. As you can
	// see, it converts from binary to decimal in reverse order. This doesn't
	// matter as long as the method returns a unique and replicable identifier.
	uint32_t _toInt(map<uint32_t, bool> history) {
		uint32_t result = 0;
		for(uint32_t i = 0; i < HIST_LEN; i++) {
			result += history[i] * pow(2, i);
		}
		return result;
	}

	// Ensures that the branch is recorded. If not, it creates the branch.
	void updateBranch(uint32_t pc) {
		bool pc_is_there = 0;
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				pc_is_there = 1;
			}
		}
		if (!pc_is_there) {
			branchList.push_back(Branch(pc));
		}
	}

	// Returns 1 if the branch exists in the vector, and 0 otherwise.
	bool exists(uint32_t pc) {
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				return 1;
			}
		}
		return 0;
	}

	// Returns the index of the branch located at the corresponding PC.
	uint32_t findBranch(uint32_t pc) {
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				return i;
			}
		}
		printf("\nError. Trying to find nonexistant branch.");
	}

	// This operates the "shift-register" for the branch history. It shifts
	// values to the right along the boolean array and adds the new value,
	// which is 1 if the branch is taken, and 0 otherwise.
	void addToHistory(uint32_t branchIndex, bool taken) {
		for (uint32_t i = HIST_LEN - 2; i > 0; i--) { // right shift
			branchList[branchIndex].history[i + 1] = branchList[branchIndex].history[i];
		}
		branchList[branchIndex].history[1] = branchList[branchIndex].history[0];
		branchList[branchIndex].history[0] = taken;
	}

	// Increments or decrements the pattern frequency counter.
	void addToPattern(uint32_t branchIndex, bool taken) {
		uint32_t patID = _toInt(branchList[branchIndex].history);
		if (taken) {
			if (branchList[branchIndex].pattern[patID] < 3) {
				branchList[branchIndex].pattern[patID]++;
			}
		} else {
			if (branchList[branchIndex].pattern[patID] > 0) {
				branchList[branchIndex].pattern[patID]--;
			}
		}
	}

	void update_execute ( 
	uint32_t pc, 
	uint32_t pc_next, 
	bool mispredict,
	bool is_brjmp, 
	uint32_t inst ) {
		uint32_t branchIndex;
		taken = (pc_next - pc) != 4; // check if branch was taken
		if ( is_brjmp ) {
			table[pc] = pc_next; // add to branch history
			updateBranch(pc); // insure that the branch is recorded
			branchIndex = findBranch(pc);
			addToPattern(branchIndex, taken); // update pattern frequency counter
			addToHistory(branchIndex, taken); // record if branch was taken
		}

		/* The following code is part of the implementation of a saturating
		   counter, which isn't used in this final version of the code, but
		   can be easily reimplemented with a change to predict_fetch. */
		if (taken) {
			if (branchList[branchIndex].counter < 3) {
				// increment saturating counter
				branchList[branchIndex].counter++;
			}
		}
		if (!taken && exists(pc)) {
			uint32_t branchIndex = findBranch(pc);
			if (branchList[branchIndex].counter > 0) {
				// decrement saturating counter
				branchList[branchIndex].counter--;
			}
		}
	}
};

  /*******************************************************/
 /****** CONSTRAINED Two-level Adaptive Predictor *******/
/************* Room for 170 total entries **************/

// HIST_LEN must be equal to 2 in order for the BP to meet the state limit
// in this particular implementation.

// Four different memory structures

class constrained_tlaBP : public BranchPredictor
{
private:

	struct Branch { // pattern frequency table
		uint32_t pc;
		std::map<uint32_t, bool> history; // taken history
		std::map<uint32_t, uint32_t> pattern; // pattern frequency table

		Branch(uint32_t _pc) {
			pc = _pc; // set pc
			// The following creates a table entry for every possible
			// pattern that is HIST_LEN bits long.
			uint32_t maxPat = 0;
			for (uint32_t i = 0; i < HIST_LEN; i++) {
				maxPat += pow(2, i);
			}
			for (uint32_t i = 0; i <= maxPat; i++) {
				pattern[i] = 0; // set initial pattern frequencies to 0
			}
		}
	};

	typedef std::vector<Branch> branches;
	branches branchList;
	std::map<uint32_t, uint32_t> table;
	bool taken;
	bool memFull;
	uint32_t index;

public:
	constrained_tlaBP ( struct bp_io& _io ) : BranchPredictor ( _io ) { }
	~constrained_tlaBP ( ) {
		taken = 0;
		memFull = 0;
		index = 0;
	}

	uint32_t newIndex(uint32_t index) {
		if (index < CONSTRAIN + 1) {
			return index++;
		} else {
			memFull = 1;
			return 0;
		}
	}

	uint32_t predict_fetch (uint32_t pc) {
		if (!exists(pc)) {
			return 0;
		}
		uint32_t prediction = predict(pc); // pattern prediction
		if (prediction == 2 || prediction == 3) { // predict taken
			return table[pc];
		}
		if (prediction == 0 || prediction == 1) { // predict not taken
			return 0;
		}
		return 0;
	}

	// Returns the two-level adaptive predictor counter.
	uint32_t predict(uint32_t pc) {
		uint32_t branchIndex = findBranch(pc);
		uint32_t patID = _toInt(branchList[branchIndex].history);
		return branchList[branchIndex].pattern[patID];
	}

	// Converts the boolean branch history into a unique integer. As you can
	// see, it converts from binary to decimal in reverse order. This doesn't
	// matter as long as the method returns a unique and replicable identifier.
	uint32_t _toInt(map<uint32_t, bool> history) {
		uint32_t result = 0;
		for(uint32_t i = 0; i < HIST_LEN; i++) {
			result += history[i] * pow(2, i);
		}
		return result;
	}

	// Ensures that the branch is recorded. If not, it creates the branch.
	void updateBranch(uint32_t pc) {
		bool pc_is_there = 0;
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				pc_is_there = 1;
			}
		}
		if (!pc_is_there) {
			if (!memFull) {
				branchList.push_back(Branch(pc));
			} else {
				branchList[index] = Branch(pc); // overwrite
			}
			index = newIndex(index);
		}
	}

	// Returns 1 if the branch exists in the vector, and 0 otherwise.
	bool exists(uint32_t pc) {
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				return 1;
			}
		}
		return 0;
	}

	// Returns the index of the branch located at the corresponding PC.
	uint32_t findBranch(uint32_t pc) {
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				return i;
			}
		}
		printf("\nError. Trying to find nonexistant branch.");
	}

	// This operates the "shift-register" for the branch history. It shifts
	// values to the right along the boolean array and adds the new value,
	// which is 1 if the branch is taken, and 0 otherwise.
	void addToHistory(uint32_t branchIndex, bool taken) {
		for (uint32_t i = HIST_LEN - 2; i > 0; i--) { // right shift
			branchList[branchIndex].history[i + 1] = branchList[branchIndex].history[i];
		}
		branchList[branchIndex].history[1] = branchList[branchIndex].history[0];
		branchList[branchIndex].history[0] = taken;
	}

	// Increments or decrements the pattern frequency counter.
	void addToPattern(uint32_t branchIndex, bool taken) {
		uint32_t patID = _toInt(branchList[branchIndex].history);
		if (taken) {
			if (branchList[branchIndex].pattern[patID] < 3) {
				branchList[branchIndex].pattern[patID]++;
			}
		} else {
			if (branchList[branchIndex].pattern[patID] > 0) {
				branchList[branchIndex].pattern[patID]--;
			}
		}
	}

	void update_execute ( 
	uint32_t pc, 
	uint32_t pc_next, 
	bool mispredict,
	bool is_brjmp, 
	uint32_t inst ) {
		uint32_t branchIndex;
		taken = (pc_next - pc) != 4; // check if branch was taken
		if ( is_brjmp ) {
			table[pc] = pc_next; // add to branch history
			updateBranch(pc); // insure that the branch is recorded
			branchIndex = findBranch(pc);
			addToPattern(branchIndex, taken); // update pattern frequency counter
			addToHistory(branchIndex, taken); // record if branch was taken
		}
	}
};

  /**********************************************/
 /****** CONSTRAINED Saturation Counter *******/
/********************************************/

// Two different (used) memory structures

class scBP : public BranchPredictor
{
private:

	struct Branch { // saturating counter and pattern frequency table
		uint32_t pc;
		uint32_t counter; // saturating counter
		std::map<uint32_t, bool> history; // taken history
		std::map<uint32_t, uint32_t> pattern; // pattern frequency table

		Branch(uint32_t _pc) {
			pc = _pc; // set pc
			counter = 0; // start counter at 0
			// The following creates a table entry for every possible
			// pattern that is HIST_LEN bits long.
			uint32_t maxPat = 0;
			for (uint32_t i = 0; i < HIST_LEN; i++) {
				maxPat += pow(2, i);
			}
			for (uint32_t i = 0; i <= maxPat; i++) {
				pattern[i] = 0; // set initial pattern frequencies to 0
			}
		}
	};
	
	typedef std::vector<Branch> branches;
	branches branchList;
	std::map<uint32_t, uint32_t> table;
	bool taken;

public:
	scBP ( struct bp_io& _io ) : BranchPredictor ( _io ) { }
	~scBP ( ) {
		taken = 0;
	}

	uint32_t predict_fetch ( uint32_t pc ) {
		if (!exists(pc)) {
			return 0;
		}
		uint32_t prediction = predict(pc);
		if (prediction == 2 || prediction == 3) { // predict taken
			return table[pc];
		}
		if (prediction == 0 || prediction == 1) { // predict not taken
			return 0;
		}
		return 0;
	}

	// Returns the two-level adaptive predictor counter.
	uint32_t predict(uint32_t pc) {
		uint32_t branchIndex = findBranch(pc);
		uint32_t patID = _toInt(branchList[branchIndex].history);
		return branchList[branchIndex].counter;
	}

	// Converts the boolean branch history into a unique integer. As you can
	// see, it converts from binary to decimal in reverse order. This doesn't
	// matter as long as the method returns a unique and replicable identifier.
	uint32_t _toInt(map<uint32_t, bool> history) {
		uint32_t result = 0;
		for(uint32_t i = 0; i < HIST_LEN; i++) {
			result += history[i] * pow(2, i);
		}
		return result;
	}

	// Ensures that the branch is recorded. If not, it creates the branch.
	void updateBranch(uint32_t pc) {
		bool pc_is_there = 0;
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				pc_is_there = 1;
			}
		}
		if (!pc_is_there) {
			branchList.push_back(Branch(pc));
		}
	}

	// Returns 1 if the branch exists in the vector, and 0 otherwise.
	bool exists(uint32_t pc) {
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				return 1;
			}
		}
		return 0;
	}

	// Returns the index of the branch located at the corresponding PC.
	uint32_t findBranch(uint32_t pc) {
		for (uint32_t i = 0; i < branchList.size(); i++) {
			if (branchList[i].pc == pc) {
				return i;
			}
		}
		printf("\nError. Trying to find nonexistant branch.");
	}

	// This operates the "shift-register" for the branch history. It shifts
	// values to the right along the boolean array and adds the new value,
	// which is 1 if the branch is taken, and 0 otherwise.
	void addToHistory(uint32_t branchIndex, bool taken) {
		for (uint32_t i = HIST_LEN - 2; i > 0; i--) { // right shift
			branchList[branchIndex].history[i + 1] = branchList[branchIndex].history[i];
		}
		branchList[branchIndex].history[1] = branchList[branchIndex].history[0];
		branchList[branchIndex].history[0] = taken;
	}

	// Increments or decrements the pattern frequency counter.
	void addToPattern(uint32_t branchIndex, bool taken) {
		uint32_t patID = _toInt(branchList[branchIndex].history);
		if (taken) {
			if (branchList[branchIndex].pattern[patID] < 3) {
				branchList[branchIndex].pattern[patID]++;
			}
		} else {
			if (branchList[branchIndex].pattern[patID] > 0) {
				branchList[branchIndex].pattern[patID]--;
			}
		}
	}

	void update_execute ( 
	uint32_t pc, 
	uint32_t pc_next, 
	bool mispredict,
	bool is_brjmp, 
	uint32_t inst ) {
		uint32_t branchIndex;
		taken = (pc_next - pc) != 4; // check if branch was taken
		if ( is_brjmp ) {
			table[pc] = pc_next; // add to branch history
			updateBranch(pc); // insure that the branch is recorded
			branchIndex = findBranch(pc);
			addToPattern(branchIndex, taken); // update pattern frequency counter
			addToHistory(branchIndex, taken); // record if branch was taken
		}

		/* The following code is part of the implementation of a saturating
		   counter, which isn't used in this final version of the code, but
		   can be easily reimplemented with a change to predict_fetch. */
		if (taken) {
			if (branchList[branchIndex].counter < 3) {
				// increment saturating counter
				branchList[branchIndex].counter++;
			}
		}
		if (!taken && exists(pc)) {
			uint32_t branchIndex = findBranch(pc);
			if (branchList[branchIndex].counter > 0) {
				// decrement saturating counter
				branchList[branchIndex].counter--;
			}
		}
	}
};

//
// Baseline Branch Predictor: simple BTB
//

#define BTB_ADDR_BITS 4
#define BTB_ENTRIES (1 << (BTB_ADDR_BITS-1))

// Sample branch predictor provided for you: a simple branch target buffer.
class BTB : public BranchPredictor 
{
private:
	typedef struct {
		uint32_t target_pc;
		uint32_t tag_pc;
	}BTBEntry_t;
	BTBEntry_t* table;

public:
	BTB ( struct bp_io& _io ) : BranchPredictor ( _io )
	{
		table = new BTBEntry_t[BTB_ENTRIES];
		memset ( table, 0, sizeof(BTBEntry_t) * BTB_ENTRIES );
	}

	~BTB ( )
	{
		delete[] table;
	}

	// Given a PC, figure out which row of the BTB table to examine
	inline uint32_t index ( const uint32_t pc )
	{
		// Extract lower BTB_ADDR_BTS bits from (pc >> 2)
		// Shift PC right two because lower two bits always zero.
		const uint32_t mask = (1 << (BTB_ADDR_BITS-1)) - 1;
		return (pc >> 2) & mask;
	}

	uint32_t predict_fetch ( uint32_t pc )
	{
		BTBEntry_t entry = table[index ( pc )];

		// Only return a prediction of the entry's tag matches this
		// PC in order to avoid aliasing.
		if ( entry.tag_pc == pc ) 
		return entry.target_pc;

		return 0;
	}

	void update_execute ( uint32_t pc, 
	uint32_t pc_next, 
	bool mispredict,
	bool is_brjmp, 
	uint32_t inst )
	{
		if( inst == BUBBLE )
		return;

		uint32_t btb_index = index ( pc );
		BTBEntry_t new_entry, 
		old_entry = table[btb_index];

		new_entry.target_pc = pc_next;
		new_entry.tag_pc = pc;

		if ( is_brjmp )
		table[btb_index] = new_entry;
	}
};

//
// A fully-associative, infinitely large BTB.  This code demonstrates how
// to use the STL map container as a fully-associative table.
//
class InfiniteBTB : public BranchPredictor
{
private:
	// Map PC to Predicted Target or zero if none (PC+4)
	std::map<uint32_t, uint32_t> table;
public:
	InfiniteBTB ( struct bp_io& _io ) : BranchPredictor ( _io ) { }
	~InfiniteBTB ( ) { }

	uint32_t predict_fetch ( uint32_t pc )
	{
		if ( table.find( pc ) != table.end() ) // Does table contain pc?
		return table[pc];  
		return 0;
	}

	void update_execute ( 
	uint32_t pc, 
	uint32_t pc_next, 
	bool mispredict,
	bool is_brjmp, 
	uint32_t inst )
	{
		if ( is_brjmp )
		table[pc] = pc_next;
	}
};

//
// No Branch Predictrion: always predict "not taken"
//
class NoBP: public BranchPredictor
{
public:
	NoBP ( struct bp_io& _io ) : BranchPredictor ( _io ) 
	{
	}

	~NoBP ( ) 
	{ 
	}

	uint32_t predict_fetch ( uint32_t pc )
	{
		return 0;
	}

	void update_execute ( 
	uint32_t pc, 
	uint32_t pc_next, 
	bool mispredict,
	bool is_brjmp, 
	uint32_t inst )
	{
	}
};




// 
// Set control signals to emulator
// (You should probably ignore this code)
//

void BranchPredictor::clock_lo ( dat_t<1> reset )
{
	if( reset.lo_word ( ) || !io.stats_reg_ptr->lo_word ( ) )
	return;

	// Examine instruction in execute stage and use it to call update_execute()
	update_execute_base (
	io.exe_reg_pc_ptr->lo_word ( ),
	io.exe_pc_next_ptr->lo_word ( ),
	io.exe_mispredict_ptr->lo_word ( ),
	// BR_N (no branch/jump) = 0 (see consts.scala)
	io.exe_br_type_ptr->lo_word ( ) != 0, 
	io.exe_reg_inst_ptr->lo_word ( ) );
}



void BranchPredictor::clock_hi ( dat_t<1> reset ) 
{
	if ( reset.lo_word ( ) ) 
	return;

	// Extract PC of instruction being fetched, and call predict_fetch with it,
	// and use the prediction to set relevant control signals back in the 
	// processor.
	uint32_t if_pc        = io.if_pc_reg_ptr->lo_word();
	uint32_t if_pred_targ = predict_fetch ( if_pc );
	*(io.if_pred_target_ptr) = LIT<32>( if_pred_targ );
	*(io.if_pred_taken_ptr) = LIT<1>( if_pred_targ != 0 );
}


void BranchPredictor::update_execute_base ( 
uint32_t pc,        // PC of this inst (in execute)
uint32_t pc_next,   // actual next PC of this inst
bool mispredict,    // Did we mispredict?
bool     is_brjmp,  // is actually a branch or jump
uint32_t inst )     // The inst itself, in case you
{                                             // want to extract arbitrary info
	if ( inst != BUBBLE ) 
	{
		brjmp_count += 1 & (long) is_brjmp;
		inst_count++;
		mispred_count += 1 & (long) ( mispredict );
	}
	cycle_count++;
	update_execute ( pc, pc_next, mispredict, is_brjmp, inst );
}


BranchPredictor::BranchPredictor ( struct bp_io& _io )  : io(_io) 
{
	brjmp_count   = 0;
	mispred_count = 0;
	inst_count    = 0;
	cycle_count   = 0;
}


BranchPredictor::~BranchPredictor ( )
{
	fprintf ( stderr, "## BRJMPs %ld\n", brjmp_count );
	fprintf ( stderr, "## INSTS %ld\n", inst_count );
	fprintf ( stderr, "## MISREDICTS %ld\n", mispred_count );
	fprintf ( stderr, "## CYCLES %ld\n", cycle_count );
}


BranchPredictor* BranchPredictor::make_branch_predictor ( struct bp_io& io )
{
	return new BRANCH_PREDICTOR ( io );
}
