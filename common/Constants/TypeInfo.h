#ifndef __TypeInfo_h
#define __TypeInfo_h

#include <string>

static const char *names[] = { "unknown", "char", "short", "unsigned_short",
                               "int", "float", "double" };
static const size_t sizes[] = { 0, 1, 2, 2, 4, 4, 8 };

class TypeInfo
{
public:
  enum TypeId
  {
    ID_UNKNOWN = 0,
    ID_CHAR,
    ID_SHORT,
    ID_USHORT,
    ID_INT,
    ID_FLOAT,
    ID_DOUBLE,
    NUM_TYPES
  };

  explicit TypeInfo(TypeId id);

  TypeId getId() const;
  const char* name() const;
  size_t size() const;

private:
  TypeId id;
};

inline TypeInfo::TypeInfo(TypeId id)
  : id(id)
{
}

inline TypeInfo::TypeId TypeInfo::getId() const
{
  return id;
}

inline const char* TypeInfo::name() const
{
  return names[id];
}

inline size_t TypeInfo::size() const
{
  return sizes[id];
}

inline TypeInfo createTypeInfo(TypeInfo::TypeId id)
{
  return TypeInfo(id);
}

inline TypeInfo createTypeInfo(const char *name)
{
  std::string nm(name);
  for (int i = 1; i < TypeInfo::NUM_TYPES; ++i)
    {
    if (nm == names[i])
      {
      return TypeInfo(static_cast<TypeInfo::TypeId>(i));
      }
    }
  return TypeInfo(TypeInfo::ID_UNKNOWN);
}

template <typename T>
TypeInfo createTypeInfo()
{
  return TypeInfo(TypeInfo::ID_UNKNOWN);
}

template<>
TypeInfo createTypeInfo<char>()
{
  return TypeInfo(TypeInfo::ID_CHAR);
}

template<>
TypeInfo createTypeInfo<short>()
{
  return TypeInfo(TypeInfo::ID_SHORT);
}

template<>
TypeInfo createTypeInfo<unsigned short>()
{
  return TypeInfo(TypeInfo::ID_SHORT);
}

template<>
TypeInfo createTypeInfo<int>()
{
  return TypeInfo(TypeInfo::ID_INT);
}

template<>
TypeInfo createTypeInfo<float>()
{
  return TypeInfo(TypeInfo::ID_FLOAT);
}

template<>
TypeInfo createTypeInfo<double>()
{
  return TypeInfo(TypeInfo::ID_DOUBLE);
}

#endif
