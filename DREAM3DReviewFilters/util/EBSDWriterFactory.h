#ifndef _ebsdwriterfactory_h
#define _ebsdwriterfactory_h

#include <map>
#include <cmath>
#include <fstream>

#include <QtCore/QDateTime>

#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/StringDataArray.h"


class EBSDWriterFactory
{
public:
  virtual ~EBSDWriterFactory();
  EBSDWriterFactory(const EBSDWriterFactory&) = delete;
  EBSDWriterFactory& operator=(const EBSDWriterFactory&) = delete;

  static EBSDWriterFactory* Instance();
  std::function<QString(QString, IDataArray::Pointer, uint64_t)> EBSDWriterFactory::getWriter(QString name);

private:
	EBSDWriterFactory();

  void initializeDataTypes();


  static EBSDWriterFactory* self;
  std::map<QString, std::function<QString(QString, IDataArray::Pointer, uint64_t)>> m_DataTypes;
};

#endif
