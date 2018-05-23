/*
 * ContextSubtract.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */

#include "ContextSubtract.h"

ContextSubtract::ContextSubtract()
:
 _fractionalSubtract(1E-9),
  _removeAll(false),
  _removeSum(false)

{
}

ContextSubtract::~ContextSubtract()
{

}


bool ContextSubtract::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
		if (strcmp(_argv[_i], "-f") == 0) {
			if (!handle_f()) return false;
		}
		if (strcmp(_argv[_i], "-A") == 0) {
			if (!handle_A()) return false;
		}
		if (strcmp(_argv[_i], "-N") == 0) {
			if (!handle_N()) return false;
		}
	}
	return ContextIntersect::parseCmdArgs(argc, argv, _skipFirstArgs);
}

bool ContextSubtract::isValidState()
{
	if (getNoHit()) {
		// Subtract will mostly use intersect's options, but -v (notHit)
		// is one of the few that aren't valid.
		_errorMsg = "\n***** ERROR: -v option is not valid for subtract. *****";
		return false;
	}

	if (getNoHit()) {
		//-r (reciprocal fraction) is also invalid.
		_errorMsg = "\n***** ERROR: -r option is not valid for subtract. *****";
		return false;
	}
	if (!ContextIntersect::isValidState()) {
		return false;
	}
	return true;
}

bool ContextSubtract::handle_f() {
    if ((_i+1) < _argc) {
    	if (isNumeric(_argv[_i+1])) {
    		_fractionalSubtract = atof(_argv[_i + 1]);
    		if (_fractionalSubtract > 0 && _fractionalSubtract <= 1.0) {
				markUsed(_i - _skipFirstArgs);
				_i++;
				markUsed(_i - _skipFirstArgs);
				return true;
    		}
    	}
    }
	_errorMsg = "\n***** ERROR: -f option must be followed by a decimal value from (0, 1.0]. *****";
	return false;
}

bool ContextSubtract::handle_A() {
   _removeAll = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextSubtract::handle_N() {
	_removeSum = true;
	markUsed(_i - _skipFirstArgs);
	return true;
}
