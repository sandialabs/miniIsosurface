#ifndef __TypeInfo_h
#define __TypeInfo_h

#include <string>

namespace util
{
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

        TypeInfo(TypeId inId)
          : id(inId)
        {}

        TypeInfo(const char * typeName)
          : id(getTypeId(typeName))
        {}

        TypeId getId() const {
            return id;
        }

        const char* name() const {
            return names[id];
        }

        size_t size() const {
            return sizes[id];
        }

        static TypeId getTypeId(const char *name) {
            std::string nm(name);
            for (int i = 1; i < TypeInfo::NUM_TYPES; ++i) {
                if (nm == names[i]) {
                    return TypeId(i);
                }
            }
            return ID_UNKNOWN;
        }

    private:
        TypeId id;
    };

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

}

#endif
