#ifndef __TypeInfo_h
#define __TypeInfo_h

#include <string>

static const char *names[] = { "unknown", "char", "short", "unsigned_short",
		"int", "float", "double" };
static const size_t sizes[] = { 0, 1, 2, 2, 4, 4, 8 };

class TypeInfo {
public:
	enum TypeId {
		ID_UNKNOWN = 0,
		ID_CHAR,
		ID_SHORT,
		ID_USHORT,
		ID_INT,
		ID_FLOAT,
		ID_DOUBLE,
		NUM_TYPES
	};

	TypeInfo(TypeId);
	TypeInfo(const char *typeName);

	TypeId getId() const;
	const char* name() const;
	size_t size() const;

private:
	TypeId id;
};

template<typename T>
TypeInfo createTemplateTypeInfo();

template<>
TypeInfo createTemplateTypeInfo<char>();

template<>
TypeInfo createTemplateTypeInfo<short>();
template<>
TypeInfo createTemplateTypeInfo<unsigned short>();

template<>
TypeInfo createTemplateTypeInfo<int>();

template<>
TypeInfo createTemplateTypeInfo<float>();

template<>
TypeInfo createTemplateTypeInfo<double>();

#endif
