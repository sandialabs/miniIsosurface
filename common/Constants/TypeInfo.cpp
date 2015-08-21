/*
 * TypeInfo.cpp
 *
 *  Created on: Aug 21, 2015
 *      Author: sjmunn
 */

#include"../includes.h"

#include"./TypeInfo.h"

TypeInfo::TypeId getTypeId(const char *name) {
	std::string nm(name);
	for (int i = 1; i < TypeInfo::NUM_TYPES; ++i) {
		if (nm == names[i]) {
			return TypeInfo::TypeId(i);
		}
	}
	return TypeInfo::ID_UNKNOWN;
}

TypeInfo::TypeInfo(TypeId inId) :
		id(inId) {
}

TypeInfo::TypeInfo(const char * typeName) : id(getTypeId(typeName)) {
}

TypeInfo::TypeId TypeInfo::getId() const {
	return id;
}

const char* TypeInfo::name() const {
	return names[id];
}

size_t TypeInfo::size() const {
	return sizes[id];
}


template<typename T>
TypeInfo createTemplateTypeInfo() {
	return TypeInfo(TypeInfo::ID_UNKNOWN);
}

template<>
TypeInfo createTemplateTypeInfo<char>() {
	return TypeInfo(TypeInfo::ID_CHAR);
}

template<>
TypeInfo createTemplateTypeInfo<short>() {
	return TypeInfo(TypeInfo::ID_SHORT);
}

template<>
TypeInfo createTemplateTypeInfo<unsigned short>() {
	return TypeInfo(TypeInfo::ID_SHORT);
}

template<>
TypeInfo createTemplateTypeInfo<int>() {
	return TypeInfo(TypeInfo::ID_INT);
}

template<>
TypeInfo createTemplateTypeInfo<float>() {
	return TypeInfo(TypeInfo::ID_FLOAT);
}

template<>
TypeInfo createTemplateTypeInfo<double>() {
	return TypeInfo(TypeInfo::ID_DOUBLE);
}
