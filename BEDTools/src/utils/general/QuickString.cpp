#include "QuickString.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "ParseTools.h"
#include "lineFileUtilities.h"

QuickString::QuickString(size_t capacity)
: _buffer(NULL),
  _currCapacity(capacity),
  _currSize(0)
{
	build();
}

QuickString::QuickString(const QuickString &qs)
:	_buffer(NULL),
	_currCapacity(qs._currCapacity),
	_currSize(0)
{
	build();
	set(qs._buffer, qs._currSize);
}

QuickString::QuickString(const char *inBuf)
{
	size_t len = strlen(inBuf);
	_currCapacity = len +1;

	build();
	set(inBuf, len);
}

QuickString::QuickString(const string &inString)
{
	size_t len = (int)inString.size();
	_currCapacity = len +1;

	build();
	set(inString.c_str(), len);
}

QuickString::QuickString(char c)
{
	_currCapacity =2;

	build();

	char buffer[2];
	buffer[0] = c;
	buffer[1] = 0;

	set(buffer, 1);
}

void QuickString::build() {
	_buffer = (char *)malloc(_currCapacity);
	clear();
}

QuickString::~QuickString(){
	free(_buffer);
}

void QuickString::clear() {
	memset(_buffer, 0, _currCapacity);
	_currSize = 0;
}

void QuickString::release() {
	free(_buffer);
	_currCapacity = DEFAULT_CAPACITY;
	build();
}

QuickString &QuickString::operator = (const char *inBuf){
	set(inBuf, strlen(inBuf));
	return *this;
}

QuickString &QuickString::operator = (const string & inBuf){
	set(inBuf.c_str(), (int)inBuf.size());
	return *this;
}

QuickString &QuickString::operator = (const QuickString & inBuf){
	set(inBuf._buffer, (int)inBuf._currSize);
	return *this;
}

QuickString &QuickString::operator = (char val) {
	clear();
	append(val);
	return *this;
}
QuickString &QuickString::operator = (int val) {
	clear();
	append(val);
	return *this;
}

QuickString &QuickString::operator = (uint32_t val) {
	clear();
	append(val);
	return *this;
}

QuickString &QuickString::operator = (size_t val) {
	clear();
	append(val);
	return *this;
}

QuickString &QuickString::operator = (float val) {
	clear();
	append(val);
	return *this;
}

QuickString &QuickString::operator = (double val) {
	clear();
	append(val);
	return *this;
}


QuickString &QuickString::operator += (const QuickString & inBuf)
{
	append(inBuf._buffer, (int)inBuf._currSize);
	return *this;
}

QuickString &QuickString::operator +=(const string &inBuf)
{
	append(inBuf.c_str(), (int)inBuf.size());
	return *this;
}

QuickString &QuickString::operator +=(char c) {

	append(c);
	return *this;
}

QuickString &QuickString::operator += (const char *inBuf)
{
	append(inBuf, strlen(inBuf));
	return *this;
}

QuickString &QuickString::operator += (int num) {
	append(num);
	return *this;
}

QuickString &QuickString::operator += (uint32_t num) {
	append(num);
	return *this;
}

QuickString &QuickString::operator += (size_t num) {
	append(num);
	return *this;
}
QuickString &QuickString::operator += (float num) {
	append(num);
	return *this;
}

QuickString &QuickString::operator += (double num) {
	append(num);
	return *this;
}

bool QuickString::operator == (const QuickString &qs) const {
	if ( _currSize != qs._currSize) {
		return false;
	}
	for (int i= _currSize-1; i > -1; i--) {
		if (_buffer[i] != qs._buffer[i]) return false;
	}
	return true;
}

bool QuickString::operator == (const string &str) const {
	if ( _currSize != str.size()) {
		return false;
	}
	for (int i= _currSize-1; i > -1; i--) {
		if (_buffer[i] != str[i]) return false;
	}
	return true;

}

bool QuickString::operator == (const char *str) const {
	size_t inLen = strlen(str);
	if (inLen != _currSize) {
		return false;
	}
	for (int i= _currSize-1; i > -1; i--) {
		if (_buffer[i] != str[i]) return false;
	}
	return true;
}


bool QuickString::operator != (const QuickString &qs) const {
	return !(*this == qs);
}

bool QuickString::operator < (const QuickString &qs) const {
	return (memcmp(_buffer, qs._buffer, max(_currSize, qs._currSize)) < 0);
}

bool QuickString::operator > (const QuickString &qs) const {
	return (memcmp(_buffer, qs._buffer, max(_currSize, qs._currSize))> 0);
}

void QuickString::set(const char *inBuf, size_t newLen) {
	reserve(newLen);
	clear();
	memcpy(_buffer, inBuf, newLen);
	_currSize = newLen;
}

void QuickString::reserve(size_t newLen) {
	newLen++; //always leave room for a null termninator.
	if (_currCapacity <= newLen) {
		while (_currCapacity <= newLen) {
			_currCapacity = _currCapacity << 1;
		}
		_buffer = (char *)realloc(_buffer, _currCapacity );
		if (_buffer == NULL) {
			fprintf(stderr, "Error: failed to reallocate string.\n");
			_currSize = 0;
			_currCapacity = 0;
			exit(1);
		}
		//initialize newly reserved memory.
		memset(_buffer + _currSize, 0, _currCapacity - _currSize);
	}
}

void QuickString::append(char c)
{
	reserve(_currSize +1);
	_buffer[_currSize] = c;
	_currSize++;
}

void QuickString::append(const char *inBuf, size_t inBufLen)
{
	reserve(_currSize + inBufLen);
	memcpy(_buffer + _currSize, inBuf, inBufLen);
	_currSize += inBufLen;
}

void QuickString::append(int num) {
	int2str(num, *this, true);
}

void QuickString::append(uint32_t num) {
	int2str((int)num, *this, true);
}

void QuickString::append(size_t num) {
	int2str((int)num, *this, true);
}

void QuickString::append(float num) {
	append(ToString(num));
}

void QuickString::append(double num) {
	append(ToString(num));
}



QuickString &QuickString::assign(const char *inBuf, size_t inBufLen)
{
	clear();
	append(inBuf, inBufLen);
	return *this;
}

void QuickString::resize(size_t newSize, char fillChar)
{
	if (newSize > _currSize) { //grow the string, pad with fillChar
		reserve(newSize);
		memset(_buffer + _currSize, fillChar, newSize -_currSize);
	} else if (newSize < _currSize) { //cut off characters from the end
		memset(_buffer + newSize, 0, _currSize - newSize);
	}
	_currSize = newSize;
}


void QuickString::substr (QuickString &newStr, size_t pos, size_t len) const
{
	if (pos >= _currSize) {
		return;
	}
	if (pos + len >= _currSize) {
		len = _currSize - pos;
	}
	newStr.set(_buffer + pos, len);
}

ostream &operator << (ostream &out, const QuickString &str) {
	out << str._buffer;
	return out;
}
