#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/SAXParser.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/sax/AttributeList.hpp>
#include <xercesc/util/XMLNetAccessor.hpp>

#include "FaceIndex.h"
#include "ShapeParser.h"
#include "Vector3.h"
#include "MathUtils.h"
#include "Dielectrics.h"

#include <vector>

void parseShapes(RWGTree *top, RWGTree *s, RWGTree *t, char *body1x3d, char *body2x3d, const Vector3 &gap, Options &options, bool (*refineEdgeFuncBody)(Edge *, void *)) {
   
  ShapeParser parser;    
  Shape body1Shape;
  Dielectric *diel1;
  parser.parse(body1x3d,body1Shape,&diel1);
  Transform tform;
  tform.setT(gap);
  body1Shape.transform(tform);
  
  Shape body2Shape;
  Dielectric *diel2;
  parser.parse(body2x3d,body2Shape,&diel2);
  
  ShapeMesh body1Mesh(body1Shape,
                      options.sharpThresh,
                      options.concaveThresh,
                      options.relativeCostAngleConcave,
                      NULL,
                      NULL);
  
  Real maxEdgeLength;
  ShapeMesh refinedBody1Mesh(body1Shape,
                             options.sharpThresh,
                             options.concaveThresh,
                             options.relativeCostAngleConcave,
                             refineEdgeFuncBody,
                             &maxEdgeLength);
  
  ShapeMesh refinedBody2Mesh(body2Shape,
                             options.sharpThresh,
                             options.concaveThresh,
                             options.relativeCostAngleConcave,
                             refineEdgeFuncBody,
                             &maxEdgeLength);
  

  ShapeMeshDielectricPairList meshDiel;
  PointEdgeVectorMap eap;
  getEdgesAtPoints(refinedBody1Mesh,eap);
  getEdgesAtPoints(refinedBody2Mesh,eap);
  meshDiel.push_back(ShapeMeshDielectricPair(&refinedBody1Mesh,diel1));
  meshDiel.push_back(ShapeMeshDielectricPair(&refinedBody2Mesh,diel2));
  Real minSegmentRadius = 0.0;
  RWGTree rwgTree(meshDiel,eap,options.minSegmentSize,minSegmentRadius);
  s = rwgTree.sub[0];
  t = rwgTree.sub[1];    
  top = &rwgTree;
}

Params attributesToParams(const AttributeList& atts) 
{
  Params params;
  for (XMLSize_t i = 0; i < atts.getLength(); i++) {
    string key = XMLString::transcode(atts.getName(i));
    if(key == "material") {
    } else {
      ParamValue val = Params::parseParam(XMLString::transcode(atts.getValue(i)));
      params.params[key] = val;
    }
  }  
  return params;
}


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
    string name = XMLString::transcode(attributes.getValue("material"));
    Params params = attributesToParams(attributes);
    dielectric = Dielectrics::getDielectric(name,params);
    if(!dielectric) abort();
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
