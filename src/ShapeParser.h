#ifndef SHAPEPARSER_H
#define SHAPEPARSER_H

#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/parsers/SAXParser.hpp>
#include "Transform.h"
#include "Vector3.h"
#include "Shape.h"
#include "Dielectric.h"
#include "RWGTree.h"
#include "Options.h"

using namespace std;
#include <vector>

XERCES_CPP_NAMESPACE_USE

class ShapeParser {
public:
  ShapeParser();
  ~ShapeParser();
  void parse(char *file, Shape &shape, Dielectric **diel);

 protected:
  SAXParser* parser;

  bool doNamespaces;
  bool doSchema;
  bool schemaFullChecking;
  SAXParser::ValSchemes valScheme;
};

class ShapeParserHandler : public HandlerBase {
public:
  ShapeParserHandler(Shape &shape);

  void startElement(const XMLCh* const name, AttributeList& attributes);
  void endElement(const XMLCh* const name);

  Dielectric *dielectric;

protected:
  Shape &shape;
  vector<Transform> tform;
};


void parseShapes(RWGTree *top, RWGTree *s, RWGTree *t, char *body1x3d, char *body2x3d, const Vector3 &gap, Options &options);

#endif
