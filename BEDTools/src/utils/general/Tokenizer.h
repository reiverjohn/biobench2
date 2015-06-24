/*
 * Tokenizer.h
 *
 *  Created on: Apr 15, 2014
 *      Author: nek3d
 */

#ifndef TOKENIZER_H_
#define TOKENIZER_H_

#include "QuickString.h"
#include <vector>

using namespace std;

class Tokenizer {
public:
	Tokenizer();
	~Tokenizer();

	// If you know the expected number of items, set this.
	// If not, don't worry about it.
	void setNumExpectedItems(int val);

	int tokenize(const QuickString &str, char delimiter = '\t', bool eofHit = false);

	// If the final element ends before a delim char, that means
	// the buffer passed in ends mid-element. The last, incomplete
	// element found can either be:
	// 1) Used now. We want it whether it's complete or not.
	// 3) Ignored altogether.
	typedef enum { USE_NOW, IGNORE } lastElemCode;
	void setKeepFinalIncompleteElem(lastElemCode code);

	//final number of valid elems may be less than total number of elems,
	//because elems are not necessarily deleted between subsequent calls
	//to tokenizer.
	int getNumValidElems() const { return _numValidElems; }
	int getNumTotalElems() const { return (int)_elems.size(); }
	const QuickString &getElem(int i) const { return (*(_elems[i])); }



private:
	static const int DEFAULT_PARSE_BUFFER_SIZE = 4096; // 8Kb
	static const int INITIAL_NUM_ELEMS = 10;
	vector<QuickString *> _elems;
	int _numExpectedElems;
	lastElemCode _keepFinalIncElem;
	int _numValidElems;

	QuickString *fetchElem(int idx);
	void resize(int newSize);
};


#endif /* TOKENIZER_H_ */
