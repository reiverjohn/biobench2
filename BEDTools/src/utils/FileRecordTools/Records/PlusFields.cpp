#include "PlusFields.h"
#include "SingleLineDelimTextFileReader.h"

PlusFields::PlusFields()
{
}

bool PlusFields::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	size_t numFields = fileReader->getNumFields() - _numOffsetFields;

	if (size() != numFields) {
		_fields.resize(numFields);
	}

	for (size_t i=0; i < numFields; i++) {
		fileReader->getField(i + _numOffsetFields, _fields[i]);
	}
	return true;
}

void PlusFields::clear() {
	//don't destroy the strings if we don't have to. Just clear their memory.
	for (int i=0; i < (int)_fields.size(); i++) {
		_fields[i].clear();
	}
}

const QuickString &PlusFields::getField(int fieldNum) const
{
	return _fields[fieldNum - _numOffsetFields - 1];
}

void PlusFields::printFields(QuickString &outBuf) const {
	for (size_t i=0; i < size(); i++) {
		outBuf.append('\t');
		outBuf.append(_fields[i]);
	}

}
