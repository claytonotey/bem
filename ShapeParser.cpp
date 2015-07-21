#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/SAXParser.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/sax/AttributeList.hpp>
#include <xercesc/util/XMLNetAccessor.hpp>

#include "defs.h"
#include "Vector3.h"
#include "FaceIndex.h"
#include "MathUtils.h"
#include "ShapeParser.h"
#include "HydrodynamicDielectric.h"
#include "Silica.h"
#include "SiliconCarbide.h"

#include <vector>
using namespace std;

ShapeParserHandler :: ShapeParserHandler(Shape &shape_) : HandlerBase(), shape(shape_)
{
}

void ShapeParserHandler :: endElement(const XMLCh* const name)
{
  tform.pop_back();
}

void ShapeParserHandler :: startElement(const XMLCh* const name, AttributeList& attributes)
{
  const XMLCh *val;
  Transform t;
  if(!XMLString::compareString(name,XMLString::transcode("Transform"))) {
    val = attributes.getValue("center");
    if(val) { Vector3 C(XMLString::transcode(val)); t.setC(C); }
    val = attributes.getValue("rotation");
    if(val) { Rotation R(XMLString::transcode(val)); t.setR(R); }
    val = attributes.getValue("scale");
    if(val) { Vector3 S(XMLString::transcode(val)); t.setS(S); }
    val = attributes.getValue("scaleOrientation");
    if(val) { Rotation SR(XMLString::transcode(val)); t.setSR(SR); }
    val = attributes.getValue("translation");
    if(val) { Vector3 T(Vector3(XMLString::transcode(val))); t.setT(T); }
  }
  tform.push_back(t);

  if(!XMLString::compareString(name,XMLString::transcode("Dielectric"))) {
    val = attributes.getValue("material");
    if(val) {
      if(!XMLString::compareString(val,XMLString::transcode("Silica")))
        dielectric = new Silica;
      if(!XMLString::compareString(val,XMLString::transcode("Au")))
        dielectric = new Au;
      if(!XMLString::compareString(val,XMLString::transcode("SiliconCarbide")))
        dielectric = new SiliconCarbide;
    }
    val = attributes.getValue("eps");
    if(val) {
      Complex eps = parseComplex(XMLString::transcode(val));      
      val = attributes.getValue("mu");
      Complex mu(1.,0.);
      if(val) {
        mu = parseComplex(XMLString::transcode(val));
      }
      dielectric = new ConstantDielectric(eps,mu);
    }
  } else if(!XMLString::compareString(name,XMLString::transcode("IndexedFaceSet"))) {
    val = attributes.getValue("coordIndex");
    if(val) {
      char *coordIndex = XMLString::transcode(val);
      int i = 0;
      while(i != -1) {
        try {
          FaceIndex index(coordIndex+i);
          shape.faceIndex.push_back(index);
        } catch(FaceIndexException) {
        }
        i = XMLString::indexOf(coordIndex,',',i);
        if(i != -1) i++;
      }      
    }
  } else if(!XMLString::compareString(name,XMLString::transcode("Coordinate"))) {
    val = attributes.getValue("point");
    if(val) {
      char *point = XMLString::transcode(val);
      int i = 0;
      while(i != -1) {
        try {
          Vector3 p(point+i);
          shape.coord.push_back(transform(tform,p));
        } catch(Vector3Exception) {
        }
        i = XMLString::indexOf(point,',',i);
        if(i != -1) i++;
      }      
    }
  }
}

ShapeParser :: ShapeParser()
{
  // Initialize the XML4C2 system
  try {
    XMLPlatformUtils::Initialize();
  }
  
  catch (const XMLException& toCatch) {
    cerr << "Error during initialization! :\n" << XMLString::transcode(toCatch.getMessage()) << "\n";
  }

  doNamespaces = false;
  doSchema = false;
  schemaFullChecking  = false;
  valScheme = SAXParser::Val_Auto;
  
  parser = new SAXParser;
  parser->setValidationScheme(valScheme);
  parser->setDoNamespaces(doNamespaces);
  parser->setDoSchema(doSchema);
  parser->setValidationSchemaFullChecking(schemaFullChecking);

}

ShapeParser :: ~ShapeParser()
{
  delete parser;
  XMLPlatformUtils::Terminate();
}

void ShapeParser :: parse(char *file, Shape &shape, Dielectric **diel)
{
  try {
    ShapeParserHandler handler(shape);
    parser->setDocumentHandler(&handler);
    parser->parse(file);
    *diel = handler.dielectric;
  }
  catch (const NetAccessorException& toCatch) {
    cerr << "NetAccessorException: " << XMLString::transcode(toCatch.getMessage()) << "\n";
    abort();
  }
  catch (const OutOfMemoryException&) {
    cerr << "OutOfMemoryException\n";
    abort();
  }
  catch (const XMLException& toCatch) {
    cerr << "An error occurred\n  Error: " << XMLString::transcode(toCatch.getMessage()) << "\n";
    abort();
  }

}

