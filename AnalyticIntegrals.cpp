#include "AnalyticIntegrals.h"

/* 

 When evaluated at a point p on the rays defined by v1 and v2 of an RWG triangle, these functions return the appropriate limiting result, so that when an identical or neighboring RWG is dotted with this result on the intersection of the domains of the two RWG, the value equals the limit as the point p approaches the boundary.

 When evaluated at a point p on the edge of an RWG, the result is in general undefined (e.g. LgradientR3).  Thus, care must be taken when evaluating these functions on the edge of an RWG, i.e. when integrating.  Note that the integration over an RWG of these functions is still well defined and non-singular, for the RWG are defined to be zero on the edge.

 */ 

int countsqrt = 0;
int countlog = 0;
int countpower = 0;

inline Real Power(Real x, Real y)
{
  return pow(x,y);
}

inline Complex Power(Complex x, Real y)
{
  return pow(x,y);
}

inline Real Sqrt(Real x)
{
  countsqrt++;
  if(x < 0)
    printf("%d sqrt < 0 %g\n",countsqrt,x);
  return sqrt(x);
}

inline Complex Sqrt(Complex x)
{
  return sqrt(x);
}

inline Real ArcTan(Real y, Real x)
{
  return atan2(y,x);
}

inline Real ArcTan(Real x)
{
  return atan(x);
}

inline Real Log(Real x)
{
  return log(x);
}

inline Complex Log(Complex x)
{
  return log(x);
}

// z = 0
void z0_S0_S1_G0_G1_G3(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real v1x, Real v1y, Real yx0, Real yy0)
{
	Real t000 = Square(x0) + Square(y0);
	Real t001 = Sqrt(t000);
	Real t002 = 1. + t000 - 2. * x0;
	Real t003 = Sqrt(t002);
	Real t004 = -1. + x0;
	Real t005 = -(y0 * Log((t003 + t004) / (t001 + x0)));
	Real t006 = Square(v1x);
	Real t007 = -1. + v1y;
	Real t008 = t006 + Square(t007);
	Real t009 = 1. / t008;
	Real t010 = -1. + yy0;
	Real t011 = Square(yx0);
	Real t012 = Square(t010) + t011;
	Real t013 = v1y * yy0;
	Real t014 = v1x * yx0;
	Real t015 = t010 - t013 - t014 + v1y;
	Real t016 = t006 + Square(v1y);
	Real t017 = t011 + Square(yy0);
	Real t018 = -2. * t013 - 2. * t014 + t016 + t017;
	Real t020 = -t013 - t014;
	Real t021 = 1. / t016;
	Real t023 = Square(v1y * yx0 - v1x * yy0);
	Real t025 = v1y - yy0;
	Real t026 = t006 - t014 + t025 * v1y;
	Real t028 = Square(t010 * v1x + yx0 - v1y * yx0);
	Real t030 = t006 - t014 + t007 * t025;
	Real t032 = Square(v2x);
	Real t033 = t032 + Square(v2y);
	Real t034 = 1. + t033 - 2. * v2x;
	Real t035 = 1. / t034;
	Real t036 = v2y * y0;
	Real t037 = v2x * x0;
	Real t038 = t004 - t036 - t037 + v2x;
	Real t039 = t000 + t033 - 2. * t036 - 2. * t037;
	Real t041 = -t036 - t037;
	Real t042 = 1. / t033;
	Real t044 = -(v2x * y0);
	Real t045 = v2y * x0;
	Real t046 = t044 + t045;
	Real t047 = Square(t046);
	Real t049 = t032 - t037 + v2y * (v2y - y0);
	Real t051 = Square(t044 + t004 * v2y + y0);
	Real t053 = t033 - t036 + x0 - v2x * (1. + x0);
	Real t055 = Sqrt(t033);
	Real t056 = 1./ (t033 * t055);
	Real t057 = Sqrt(t017);
	Real t058 = Sqrt(t016);
	Real t059 = 1. / (t016 * t058);
	Real t060 = 1. / t055;
	Real t061 = Sqrt(t012);
	Real t062 = Sqrt(t008);
	Real t063 = 1. / (t008 * t062);
	Real t064 = Sqrt(t034);
	Real t065 = 1. / (t034 * t064);
	Real t066 = 1. / t064;
	Real t067 = Sqrt(t039);
	Real t068 = Sqrt(t018);
	Real t069 = Log(-((t033 - t036 - t037 + t055 * t067) / (-t041 - t001 * t055)));
	Real t070 = Log(-((-t013 - t014 + t016 + t058 * t068) / (-t020 - t057 * t058)));
	Real t071 = Log((t053 + t064 * t067) / (t038 + t003 * t064));
	Real t072 = Log((-t014 + t016 + t062 * t068 + yy0 - v1y * (1. + yy0)) / (t015 + t061 * t062));
  
  * S0 = 0.5;
  * S1 = (t005 + t046 * t060 * t069 + t066 * t071 * (-t045 + v2y + (-1. + v2x) * y0)) / v2y;
  * G0x = (-1. - v2x + 3. * x0) / 6.;
  * G0y = (-1. - v1y + 3. * yy0) / 6.;
  * G1x = (t003 * t035 * t038 - t001 * t041 * t042 + t042 * t049 * t067 - t035 * t053 * t067 + t047 * t056 * t069 - t051 * t065 * t071) / 2.;
  * G1y = (-(t020 * t021 * t057) + t009 * t015 * t061 + t021 * t026 * t068 - t009 * t030 * t068 + t023 * t059 * t070 - t028 * t063 * t072) / 2.;
  * G3x = -(t060 * t069) + t066 * t071;
  * G3y = -(t070 / t058) + t072 / t062;  
}

// z = 0
void z0_S0_S1_G0_G1(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real v2x, Real v2y, Real x0, Real y0, Real v1x, Real v1y, Real yx0, Real yy0)
{
  countlog = 0;
	countsqrt = 0;
  //printf("%g %g %g %g %g %g %g %g\n",v2x,v2y,x0,y0,v1x,v1y,yx0,yy0);
  Real t000 = Square(x0) + Square(y0);
	Real t001 = Sqrt(t000);
	Real t002 = 1. + t000 - 2. * x0;
	Real t003 = Sqrt(t002);
	Real t004 = -1. + x0;
	Real t005 = -(y0 * Log((t003 + t004) / (t001 + x0)));
	Real t006 = Square(v1x);
	Real t007 = -1. + v1y;
	Real t008 = t006 + Square(t007);
	Real t009 = 1. / t008;
	Real t010 = -1. + yy0;
	Real t011 = Square(yx0);
	Real t012 = Square(t010) + t011;
	Real t013 = v1y * yy0;
	Real t014 = v1x * yx0;
	Real t015 = t010 - t013 - t014 + v1y;
	Real t016 = t006 + Square(v1y);
	Real t017 = t011 + Square(yy0);
	Real t018 = -2. * t013 - 2. * t014 + t016 + t017;
	Real t020 = -t013 - t014;
	Real t021 = 1. / t016;
	Real t023 = Square(v1y * yx0 - v1x * yy0);
	Real t025 = v1y - yy0;
	Real t026 = t006 - t014 + t025 * v1y;
	Real t028 = Square(t010 * v1x + yx0 - v1y * yx0);
	Real t030 = t006 - t014 + t007 * t025;
	Real t032 = Square(v2x);
	Real t033 = t032 + Square(v2y);
	Real t034 = 1. + t033 - 2. * v2x;
	Real t035 = 1. / t034;
	Real t036 = v2y * y0;
	Real t037 = v2x * x0;
	Real t038 = t004 - t036 - t037 + v2x;
	Real t039 = t000 + t033 - 2. * t036 - 2. * t037;
	Real t041 = -t036 - t037;
	Real t042 = 1. / t033;
	Real t044 = -(v2x * y0);
	Real t045 = v2y * x0;
	Real t046 = t044 + t045;
	Real t047 = Square(t046);
	Real t049 = t032 - t037 + v2y * (v2y - y0);
	Real t051 = Square(t044 + t004 * v2y + y0);
	Real t053 = t033 - t036 + x0 - v2x * (1. + x0);
	Real t055 = Sqrt(t033);
	Real t056 = 1./ (t033 * t055);
	Real t057 = Sqrt(t017);
	Real t058 = Sqrt(t016);
	Real t059 = 1. / (t016 * t058);
	Real t060 = 1. / t055;
	Real t061 = Sqrt(t012);
	Real t062 = Sqrt(t008);
	Real t063 = 1. / (t008 * t062);
	Real t064 = Sqrt(t034);
	Real t065 = 1. / (t034 * t064);
	Real t066 = 1. / t064;
	Real t067 = Sqrt(t039);
	Real t068 = Sqrt(t018);
	Real t069 = Log(-((t033 - t036 - t037 + t055 * t067) / (-t041 - t001 * t055)));
	Real t070 = Log(-((-t013 - t014 + t016 + t058 * t068) / (-t020 - t057 * t058)));
	Real t071 = Log((t053 + t064 * t067) / (t038 + t003 * t064));
	Real t072 = Log((-t014 + t016 + t062 * t068 + yy0 - v1y * (1. + yy0)) / (t015 + t061 * t062));
  
  * S0 = 0.5;
  * S1 = (t005 + t046 * t060 * t069 + t066 * t071 * (-t045 + v2y + (-1. + v2x) * y0)) / v2y;
  * G0x = (-1. - v2x + 3. * x0) * ONESIXTH;
  * G0y = (-1. - v1y + 3. * yy0) * ONESIXTH;
  * G1x = (t003 * t035 * t038 - t001 * t041 * t042 + t042 * t049 * t067 - t035 * t053 * t067 + t047 * t056 * t069 - t051 * t065 * t071) * 0.5;
  * G1y = (-(t020 * t021 * t057) + t009 * t015 * t061 + t021 * t026 * t068 - t009 * t030 * t068 + t023 * t059 * t070 - t028 * t063 * t072) * 0.5;
}

// z = 0

Real z0_gdotG3(Real v2x, Real v2y)
{
	Real u000 = Square(v2y);
	Real u001 = Square(v2x);
	Real u002 = -2. * v2x;
	Real u003 = u000 + u001;
	Real u004 = 1. + u002 + u003;
	Real u005 = Sqrt(u004);
	Real u006 = -1. + v2x;
	Real u007 = u005 + u006;
	Real u008 = -Log(u007 / u000);
	Real u009 = Square(u003);
	Real u010 = Square(u000);
	Real u011 = Square(u001);
	Real u012 = v2x * u001;
	Real u013 = -2. * u012;
	Real u014 = 2. * u001;
	Real u015 = u000 * u002;
	Real u016 = Sqrt(u003 + u010 + u011 + u013 + u000 * u014 + u015);
	Real u017 = -6. * u009 * u016 * Log(u000);
	Real u018 = Sqrt(u003);
	Real u019 = u018 - v2x;
	Real u020 = 1. + u005 - v2x;
	Real u021 = 6. * Log((u007 * u019) / u020);
	Real u022 = Square(u012);
	Real u023 = u012 * u001;
	Real u024 = 1. / u009;
	Real u025 = 1. / u003;
	Real u026 = 1. / u018;
	Real u027 = 1. / (u018 * u009);
	Real u028 = 1. / u004;
	Real u029 = Square(u028);
	Real u030 = 1. / u005;
	Real u031 = u029 / u005;
	Real u032 = u029 * u005;
	Real u033 = 1. / u016;
	Real u034 = 3. * v2x;
	Real u035 = 6. * u001;
	Real u036 = -3. * u000;
	Real u037 = 1. + u000;
	Real u038 = 2. * u000;
	Real u039 = 5. * u000;
	Real u040 = Log(u018 + v2x);
	Real u041 = Log(u007);
	Real u042 = Log(u019);
	Real u043 = 4. + u003 - 5. * v2x;
	Real u044 = u005 * u018;
	Real u045 = u003 - v2x;
	Real u046 = u044 + u045;
	Real u047 = Log(u046);
	Real u048 = Log(2. * u004);
	Real u049 = Log(-u000 - u001 + u044 + v2x);
	Real u050 = Log(u046 / u019);
	Real u051 = Log(u046 / u007);
	Real u052 = 1. / (u016 + u045);
	Real u053 = Log((-u001 + u000 * u006 + u012 + u016 * u018) * u025 * u052);
	Real u054 = u005 * u016;
	Real u055 = Log(u028 * u052 * (-u012 + u014 + u054 - u037 * v2x));
	Real u056 = 3. * u005 * u024 * (u003 + u034) + u027 * (u035 + u036) * u050;
	Real u057 = u001 - u038;
	Real u058 = -3. * u018 * u029;
	Real u059 = -9. * u006 * u029;
	Real u060 = -2. * u018 * u028;
	Real u061 = -9. * u024 * v2x;
	Real u062 = 3. * u031 * u051;
	Real u063 = 9. * u000 * u016;
	Real u064 = -5. - 6. * u005 - 6. * u018 + u021 - 9. * u030 - 18. * (-1. + u042) + u026 * (11. - 6. * u047) - 18. * u026 * (-1. + u047) + 9. * u026 * (-3. + 2. * u047) + u030 * (11. - 6. * u048) + u030 * (-2. + 6. * u048) - 18. * u030 * (-1. + u049) + u030 * (-11. + 6. * u049) + 18. * u026 * (-1. + u042 + u003 * u033 * u053) + 18. * u030 * (-1. + u008 - u004 * u033 * u055) + u056 + 9. * (1. + 2. * u028 - 2. * u040 + 2. * u006 * u032 * u051 + u060) + u061 + 2. * (-2. + 6. * u040 - u043 * u058 - u059 + u062 * (2. - u000 + u014 - 4. * v2x)) + u026 * (-11. + 6. * u042 + 3. * u029 * u033 * (u001 * (6. + 25. * u000 + 6. * u010) + (20. + 6. * u000) * u011 + (9. + 10. * u000) * u013 + 2. * u022 - 10. * u023 + u015 * (9. + u039) + u000 * (6. + 2. * u010 + u039)) * u053 - u058 * (5. + u014 + u038 - 7. * v2x) - 9. * u003 * u029 * (2. + u003 - 3. * v2x)) + 2. * (-1. + 3. * u040) * v2x + 2. * (-1. + 3. * u041) * v2x + (2. - 6. * u041 + u056 + u061) * v2x + v2x * (2. - 6. * u040 + u043 * u058 + u059 + u062 * (-2. + u000 - 2. * u001 + 4. * v2x)) + 3. * (u000 * u027 * u034 * u050 + u024 * u057 + u005 * (2. - u001 * u024 + u024 * u038 - u025 * v2x)) + 9. * u026 * (3. + 2. * u003 * u028 - 2. * u042 + u060 - 2. * u028 * u033 * u053 * (u000 * (2. + u000) + u011 - 3. * u012 + u014 * u037 + u036 * v2x)) + 3. * (u006 * u031 * u036 * u051 + u029 * (1. + u002 + u057) + u018 * (2. + u006 * u028 + u029 * (-1. - u001 + u038 + 2. * v2x))) + u024 * u030 * u033 * (9. * u001 * u016 + 4. * u000 * u001 * u016 + 2. * u010 * u016 + 2. * u011 * u016 + 9. * u012 * u016 + u017 + 12. * u000 * u001 * u016 * u041 + 6. * u010 * u016 * u041 + 6. * u011 * u016 * u041 + 6. * u000 * u054 + u035 * u054 + 3. * u000 * u001 * u055 + 3. * u010 * u055 + 18. * u001 * u010 * u055 + 18. * u000 * u011 * u055 - 6. * u012 * u055 - 12. * u000 * u012 * u055 + 6. * u022 * u055 - 6. * u023 * u055 + u035 * u055 + u036 * u055 - u063 - 9. * u016 * v2x + 9. * u054 * v2x + 12. * u000 * u055 * v2x - 6. * u010 * u055 * v2x + u063 * v2x + 6. * u055 * u010 * u000) + 18. * (-1. + Log(u020));

  return u064 / 36.;
}

inline Complex lognegconj(Complex &z1)
{
  Real i = imag(z1);
  if(i<0.)
    i = -PI - i;
  else
    i = PI - i;
  return Complex(real(z1),i);
}

// z != 0
void S0_S1_G0_G1_C0_C1(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real *C0z, Real *C1z, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0)
{
  countlog = 0;
	countsqrt = 0;
  //printf("%g %g %g %g %g %g %g %g %g %g\n",v2x,v2y,x0,y0,z0,v1x,v1y,yx0,yy0,yz0);
	Real t000 = -1. + x0;
	Real t001 = v2x * x0;
	Real t002 = t000 - t001 + v2x;
	Real t003 = v2y * y0;
	Real t004 = t002 - t003;
	Real t005 = 2. * t004;
	Real t006 = -1. + v2x;
	Real t007 = t001 + t003;
	Real t008 = -2. * t007;
	Real t009 = -2. * t003;
	Real t010 = -2. * x0;
	Real t011 = t010 * v2x;
	Real t012 = -2. * y0;
	Real t014 = -1. + v1y;
	Real t015 = -2. * yy0;
	Real t016 = -1. + yy0;
	Real t017 = v1x * yx0;
	Real t018 = -2. * t017;
	Real t019 = t015 * v1y;
	Real t020 = 1. + t017 + t016 * v1y - yy0;
	Real t021 = v1y * yy0;
	Real t022 = t017 + t021;
	Real t023 = v1y - yy0;
	Real t024 = 1. + x0;
	Real t025 = v2y - y0;
  Real t030 = v2x * v2y;
	Real t031 = t008 * v2x;
	Real t032 = t001 * v2y;
	Real t033 = -8. * t001;
	Real t034 = t005 * t006;
	Real t035 = 2. * t002;
	Real t036 = 8. * t002;
	Real t037 = Square(t005*t005) * t006;
	Real t038 = Square(t008*t008) * v2x;
	Real t039 = 16. * t002;
	Real t040 = 6. * t002;
	Real t041 = -16. * t001;
	Real t042 = -6. * t001;
	Real t043 = Square(z0);
	Real t044 = Square(y0);
	Real t045 = t006 * t044;
	Real t046 = -(t043 * v2x);
	Real t047 = t043 + t045 + t046;
	Real t048 = Square(v2y);
	Real t049 = -(t043 * t048);
	Real t050 = Square(t006);
	Real t051 = v2y * t048;
	Real t052 = Square(t000);
	Real t053 = -(t002 * t012 * t048) + t051 * t052;
	Real t054 = Square(t009);
	Real t055 = t043 + t044;
	Real t056 = -4. * t048 * t055;
	Real t057 = t054 + t056;
	Real t058 = t044 * v2x;
	Real t059 = t046 + t058;
	Real t060 = Square(v2x);
	Real t061 = Square(x0);
	Real t062 = t051 * t061 + t011 * t048 * y0;
	Real t064 = Square(t048);
	Real t065 = t055 + t061;
	Real t066 = t048 + t060;
	Real t067 = 2. * t066;
	Real t068 = 1. + t055 + (-2. + x0) * x0;
	Real t069 = t048 + t050;
	Real t070 = 2. * t069;	
	Real t072 = 1. / t043;
	Real t073 = t043 + t061;
	Real t076 = Sqrt(t069);
	Real t077 = 1. / t076;
	Real t083 = Square(v1x);
	Real t084 = Square(v1y);
	Real t085 = t083 + t084;
	Real t086 = Square(yy0);
	Real t087 = Square(yz0);
	Real t088 = t087 + Square(yx0);
	Real t089 = t086 + t088;
	Real t090 = t018 + t019 + t085 + t089;
	Real t091 = Sqrt(t090);
	Real t092 = 1. + t015 + t089;
	Real t093 = Sqrt(t092);
	Real t094 = Sqrt(t089);
	Real t095 = Square(t014);
	Real t096 = t083 + t095;
	Real t098 = Sqrt(t096);
	Real t097 = 1. / (t096 * t098);
	Real t100 = Sqrt(t085);
	Real t099 = 1. / (t085 * t100);
	Real t114 = Sqrt(t066);
	Real t106 = 1. / (t066 * t114);
	Real t107 = 1. / (t069 * t076);
	Real t108 = 1. / t051;
	Real t109 = 1. / t064;
	Real t110 = t006 * t050;
	Real t141 = Square(t005);
	Real t111 = t005 * t141;
	Real t113 = 1. / t048;
	Real t115 = 1. / t114;
	Real t144 = Square(t008);
	Real t116 = t008 * t144;
	Real t117 = v2x * t060;
	Real t118 = Sqrt(t043 * t048);
	Real t119 = 1. / t118;
	Real t120 = t008 + t065 + t066;
	Real t121 = Sqrt(t120);
	Real t122 = t005 + t068 + t069;
	Real t123 = Sqrt(t122);
	Real t124 = Sqrt(t065);
	Real t125 = Sqrt(t068);
	Real t127 = t048 * t073;
	Real t128 = t003 * t011 + t055 * t060 + t127;
	Real t130 = -(t002 * t009) + t050 * t055 + t048 * (1. + t010 + t073);
	Real t132 = 1. / t066;
	Real t133 = 1. / t069;
	Real t134 = t086 + t087;
	Real t135 = t014 * t016 * t018 + t088 * t095 + t083 * (1. + t015 + t134);
	Real t137 = t017 * t019 + t084 * t088 + t083 * t134;
	Real t139 = 1. / t085;
	Real t140 = 1. / t096;
	Real t142 = 4. * t068 * t069;
	Real t145 = 4. * t065 * t066;
	Real t149 = t008 + t067;
	Real t151 = t005 + t070;
	Real t153 = -t003 + t066;
	Real t154 = t153 - t024 * v2x + x0;
	Real t156 = -t001 + t060 + t025 * v2y;
  Real t168 = 2. * t114 * t124;
	Real t169 = t043 * t051 * y0;
	Real t170 = t043 * v2y;
	Real t171 = t048 * t065;
	Real t172 = 2. * t171;
	Real t173 = -(t066 * x0);
	Real t174 = -(t067 * x0);
	Real t175 = 4. * t117;
	Real t176 = -(t175 * x0);
	Real t177 = t033 * t066 + t176;
	Real t178 = t066 * v2x;
	Real t179 = -4. * t178;
	Real t180 = -4. * t061 * t066;
	Real t182 = t116 * v2x;
	Real t183 = ArcTan(t003 * t119);
	Real t184 = 4. * t118 * t183;
	Real t185 = t031 - 2. * t173;
	Real t186 = t076 * t123;
	Real t187 = t060 * v2y;
	Real t188 = t043 * t064;
	Real t190 = t043 * x0;
	Real t191 = t048 * t068;
	Real t192 = -2. * t191;
	Real t193 = -(t000 * t069);
	Real t194 = t036 * t069;
	Real t195 = -(t000 * t070);
	Real t196 = 4. * t052 * t069;
	Real t197 = t006 * t111;
	Real t198 = 4. * t110;
	Real t199 = -(t000 * t198);
	Real t200 = t006 * t069;
	Real t201 = -4. * t200;
	Real t203 = 2. * t006 * t064 * Square(t068);
	Real t204 = 2. * t064 * Square(t065) * v2x;
	Real t205 = t076 * t125;
	Real t206 = t114 * t121;
	Real t207 = t003 - t066;
	Real t208 = -(t048 * t068);
	Real t209 = t050 + t069;
	Real t210 = t141 * (-t110 + 3. * t200) - t201 * (t052 * t070 + t208) - 4. * t005 * t193 * t209;
	Real t211 = -(t048 * t065);
	Real t212 = t060 + t066;
	Real t213 = t144 * (-t117 + 3. * t178) - t179 * (t061 * t067 + t211) - 4. * t008 * t173 * t212;
	Real t214 = -(t000 * t208);
	Real t215 = 4. * t193;
	Real t216 = -(t000 * t052 * t069);
	Real t217 = -4. * t000 * Square(t069);
	Real t218 = 11. * t002 * t069 - 13. * t000 * t110;
	Real t219 = 4. * t173;
	Real t220 = -(t066 * x0 * t061);
	Real t221 = -13. * t117 * x0;
	Real t222 = -11. * t001 * t066 + t221;
	Real t223 = 4. * Square(t066) * x0;
	Real t224 = t141 - t142;
	Real t225 = t144 - t145;
	Real t226 = -ArcTan(t025 * t119 * v2y);
	Real t227 = -4. * t113 * t118 * t226;
	Real t228 = t034 - 2. * t193;
	Real t229 = t009 * t113;
	Real t230 = 4. * t009 * t109 * t118 * t226;
	Real t231 = t009 * t184;
	Real t232 = 2. * t229;
	Real t234 = t020 * t140;
	Real t235 = -t017 + t083;
	Real t236 = t235 + t023 * v1y;
	Real t238 = t014 * t023 + t235;

  // negative powers
	Complex t240 = Sqrt(Complex(t049));
	Complex t112 = t049 * t240;
	Complex t241 = Sqrt(Complex(t057));
	Complex t189 = t072 * t112 * x0;
	Complex t242 = t240 * y0;
	Complex t243 = -2. * t050 * t242;
	Complex t246 = -2. * t060 * t242;
	Complex t302 = -6. * t240;
	Complex t303 = -2. * t066 * t240;
	Complex t306 = -t001 + t241;
	Complex t307 = -t182 + t144 * (t173 + (-5. * t001 + t241) * v2x) - t172 * (t174 + (t011 + t241) * v2x) - t180 * (t173 + t306 * v2x) + t008 * (4. * t171 * v2x + (t177 - t060 * t302 - t303) * x0);
	Complex t308 = 1. / t240;
	Complex t247 = Sqrt(v2y * (t062 + t246 + t030 * (t059 + t241 * x0)));
	Complex t310 = conj(t247);
	Complex t312 = t127 * t240;
	Complex t313 = t169 + t312;
	Complex t314 = 1. / t310;
	Complex t315 = -(t308 * t314) * 0.5;
	Complex t317 = -(t121 * t240);
	Complex t318 = t024 * t242;
	Complex t320 = t189 + t048 * t061 * t240;
	Complex t321 = t320 + t003 * t240 * x0;
	Complex t322 = t002 + t241;
	Complex t323 = t070 * t240;
	Complex t324 = -t197 + t141 * (t193 + t006 * (5. * t002 + t241)) + t192 * (t195 + t006 * (t035 + t241)) + t196 * (t193 + t006 * t322) + t005 * (4. * t006 * t191 + t000 * (t194 + t199 - t050 * t302 + t323));
	Complex t328 = -(t125 * t240);
	Complex t329 = t052 * t240;
	Complex t330 = t002 * t003 * t240;
	Complex t331 = t043 + t240;
	Complex t332 = 14. * t240;
	Complex t333 = 10. * t240;
	Complex t334 = -t037 + t111 * (t193 + t006 * (7. * t002 + t241)) - 2. * (t203 + t216 * (t069 * (t040 + t241) + t050 * (t035 - t302)) + t214 * (t069 * (t036 + t241) + 2. * t050 * t322)) + t141 * (5. * t006 * t191 + t000 * (t218 + t323 + t050 * t333)) + t005 * (t208 * (t215 + t006 * (t039 - t302)) + t052 * (t217 + t200 * (24. * t002 + t332) + t110 * (4. * t002 + t333)));
	Complex t337 = -t038 + t116 * (t173 + (-7. * t001 + t241) * v2x) + t008 * (t061 * (-t223 + t178 * (-24. * t001 + t332) + t117 * (-4. * t001 + t333)) + t211 * (t219 + (t041 - t302) * v2x)) - 2. * (t204 + t220 * (t066 * (t042 + t241) + t060 * (t011 - t302)) - t211 * (t066 * (t033 + t241) + 2. * t060 * t306) * x0) + t144 * (5. * t171 * v2x + (t222 - t303 + t060 * t333) * x0);
	Complex t244 = Sqrt(v2y * (t053 + t243 + t006 * (t047 + t000 * t241) * v2y));
	Complex t338 = conj(t244);
	Complex t339 = 1. / t338;
	Complex t340 = t187 * (-t170 + t242);
	Complex t341 = (t113 * t308) * 0.5;

  // logs
  Complex z1n = (-4. * (-t188 + t313 + t310 * t317 + t340 + t032 * (t170 - t240 * (v2y + y0))));
  Complex z1d = (t060 + t207 + t240);
	Complex z1 = Log(z1n/z1d);
  Complex z2 = lognegconj(z1);

  Complex z3n = -4. * (t169 + t048 * (t043 * (-t002 + t240) + t329) + t330 + t328 * t338);
  Complex z3d = (t003 + t240);
	Complex z3 = Log(z3n/z3d);
  Complex z4 = lognegconj(z3);

  Complex z5n = (-4. * (-t112 + t169 - t188 - t048 * t190 + t321 + t317 * t338 + t340 + t030 * (-t318 + v2y * (t190 + t331 - t240 * x0))));
  Complex z5d = (t003 + t050 - t069 + t240);
	Complex z5 = Log(z5n/z5d);
  Complex z6 = lognegconj(z5);

  Complex z7n = -4. * (t032 * (t170 - t242) - t124 * t240 * t310 + t313);
  Complex z7d = (t003 + t240);
	Complex z7 = Log(z7n/z7d);
  Complex z8 = lognegconj(z7);

  Complex z1c = t307 * t315;
  Complex z2c = conj(z1c);
  Complex z3c = 0.5 * t308 * t324 * t339;
  Complex z4c = conj(z3c);
  Complex z5c = z3c;
  Complex z6c = z4c;
  Complex z7c = z1c;
  Complex z8c = z2c;

	Real t362 = Log((-t017 + t085 + t091 * t098 + yy0 - v1y * (1. + yy0)) / (t016 - t017 - t021 + t093 * t098 + v1y));
	Real t363 = Log((t022 - t085 - t091 * t100) / (t022 - t094 * t100));
	Real t367 = Log((t076 * t121 + t154) / (t004 + t205));
	Real t368 = Log((t001 - t206 + t207) / (t007 - t114 * t124));
	Real t370 = Log(t008 + 2. * (t066 + t206));
	Real t371 = Log(t005 + 2. * t205);
	Real t372 = Log(t121 + v2x - x0);
	Real t374 = Log(t005 + 2. * (t069 + t186));
	Real t375 = Log(-t000 + t006 + t123);
	Real t376 = Log(t008 + t168);
	Real t380 = 2. * t115 * t376;
  
  Real t369 = Log(t008 - t011 + t048 + t055);
	Real t373 = Log(t005 - t035 + t048 + t055);
	Real t377 = Log(t055);

  Real t383_1 = t184 + (-t008 + t011) * t377 + t185 * t380;
	Real t383 = t383_1 + real(z8c * z8) + real(z7c * z7);
	Real t384 = 4. * t375;
	Real t385_1 = (-4. + t227 +  t229 * t373 - 2. * t077 * t113 * t228 * t374 + t384) * t048;
  Real t385 = t385_1 + real(z5c * z5) + real(z6c * z6);
	Real t386 = 4. * t372;
  Real t387 = (t107 * (-2. * t005 * t205 + t224 * t371) * y0) * 0.125;
  Real t388 = (t106 * (t008 * t168 - t225 * t376) * y0) * 0.125;

  Real z11 = 0.25 * (-v2y * t109 * real(t315 * t337 * z1) + 2. * y0 * real(t307 * t314 * t341 * z1));
  Real t392 = z11 + (-(v2y * (-2. - 2. * t113 * t121 * t132 * t185 + t230 + t232 + t106 * t109 * t213 * t370 + t386 - t109 * t369 * (2. * t060 * t065 - 2. * t065 * t066 + t144 + 2. * t061 * t212 + 4. * t031 * x0))) * 0.125 + ((-4. + t227 + t229 * t369 - 2. * t113 * t115 * t185 * t370 + t386) * y0) * 0.25);

  Real z31 = 0.125 * (t108 * real(t308 * t334 * t339 * z3) + 2. * y0 * t113 * real(t308 * t324 * t339 * z3));
  Real t389 = z31 - (t108 * (-2. * t048 * t125 * t133 * t228 + t231 + t107 * t210 * t371 + (-4. * t000 * t034 - 2. * t050 * t068 + t068 * t070 - t141 - 2. * t052 * t209) * t377)) * 0.125 - (t113 * (t184 + 2. * t077 * t228 * t371 + (-t005 + t035) * t377) * y0) * 0.25;
  
  Real z51 = 0.125 * (t109 * v2y * real(t308 * t334 * t339 * z5) + 4. * y0 * t113 * real(z5c * z5));
  Real t391 = z51 - ((-2. - 2. * t113 * t123 * t133 * t228 + t230 + t232 - t109 * (4. * t000 * t034 + 2. * t050 * t068 - 2. * t068 * t069 + t141 + 2. * t052 * t209) * t373 + t107 * t109 * t210 * t374 + t384) * v2y) * 0.125 + (t113 * t385_1 * y0) * 0.25;

  Real z71 = -0.25 * (t108 * real(t337 * t315 * z7) + 2. * y0 * t113 * real(z7c * z7));
  Real t390 = z71 - (t108 * (-2. * t048 * t124 * t132 * t185 + t231 + t106 * t213 * t376 + t377 * (-2. * t060 * t065 + t065 * t067 - t144 - 2. * t061 * t212 - 4. * t031 * x0))) * 0.125 - (t113 * t383_1 * y0) * 0.25;

  Real t393 = (t107 * y0 * (2. * t151 * t186 - t224 * Log(t151 + 2. * t186))) * 0.125;
  Real t394 = (t106 * y0 * (2. * t149 * t206 - t225 * Log(t149 + 2. * t206))) * 0.125;

  // S1 correction
  Real S1corr;
  if((imag(z1c) >= 0.) == (imag(z5c) >= 0.))
    S1corr = ((imag(z1) + imag(z2) - imag(z7) - imag(z8)) + (imag(z5) + imag(z6) - imag(z3) - imag(z4)));
  else
    S1corr = ((imag(z1) + imag(z2) - imag(z7) - imag(z8)) - (imag(z5) + imag(z6) - imag(z3) - imag(z4)));
  S1corr *= imag(z1c);

  * S0 = 0.5;
  * S1 = t113 * S1corr * 0.25 + (4. + 4. * t113 * t118 * t226 + t113 * (real(z1c * z1) + real(z2c * z2)) - t229 * t369 + 2. * t113 * t115 * t185 * t370 - 4. * t372) * 0.25 - (t113 * (-4. * t118 * t183 + real(z3c * z3) + real(z4c * z4) + 2. * t077 * (-t034 + t195) * t371 + t009 * t377)) * 0.25 - (t113 * t383) * 0.25 + (t113 * t385) * 0.25;
  * G0x = (-1. - v2x + 3. * x0) * ONESIXTH;
  * G0y = (-1. - v1y + 3. * yy0) * ONESIXTH;
  * G1x = (t007 * t124 * t132 + t004 * t125 * t133 - t121 * t133 * t154 + t121 * t132 * t156 - t107 * t130 * t367 + t106 * t128 * t368) * 0.5;
  * G1y = (t022 * t094 * t139 - t093 * t234 + t091 * t139 * t236 - t091 * t140 * t238 - t097 * t135 * t362 + t099 * t137 * t363) * 0.5;
  * C0z = (-(v2y * x0) + y0 + v2x * y0) * ONESIXTH;
  * C1z = t387 + t388 + t393 - t394 + x0 * (-t389 + t390 + t391 - t392);
}

// z != 0
void S0_S1_G0_G1(Real *S0, Real *S1, Real *G0x, Real *G0y, Real *G1x, Real *G1y, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0)
{
  //printf("%g %g %g %g %g %g %g %g\n",v2x,v2y,x0,y0,v1x,v1y,yx0,yy0);
	Real t000 = -1. + x0;
	Real t001 = v2x * x0;
	Real t002 = t000 - t001 + v2x;
	Real t003 = v2y * y0;
	Real t004 = t002 - t003;
	Real t005 = 2. * t004;
	Real t006 = -1. + v2x;
	Real t007 = t001 + t003;
	Real t008 = -2. * t007;
	Real t009 = -2. * t003;
	Real t010 = -2. * x0;
	Real t011 = t010 * v2x;
	Real t012 = -2. * y0;
	Real t014 = -1. + v1y;
	Real t015 = -2. * yy0;
	Real t016 = -1. + yy0;
	Real t017 = v1x * yx0;
	Real t018 = -2. * t017;
	Real t019 = t015 * v1y;
	Real t020 = 1. + t017 + t016 * v1y - yy0;
	Real t021 = v1y * yy0;
	Real t022 = t017 + t021;
	Real t023 = v1y - yy0;
	Real t024 = 1. + x0;
	Real t025 = v2y - y0;
	Real t030 = v2x * v2y;
	Real t031 = t008 * v2x;
	Real t032 = t001 * v2y;
	Real t033 = -8. * t001;
	Real t034 = t005 * t006;
	Real t035 = 2. * t002;
	Real t036 = 8. * t002;
	Real t043 = Square(z0);
	Real t044 = Square(y0);
	Real t045 = t006 * t044;
	Real t046 = -(t043 * v2x);
	Real t047 = t043 + t045 + t046;
	Real t048 = Square(v2y);
	Real t049 = -(t043 * t048);
	Real t050 = Square(t006);
	Real t051 = v2y * t048;
	Real t052 = Square(t000);
	Real t053 = -(t002 * t012 * t048) + t051 * t052;
	Real t054 = Square(t009);
	Real t055 = t043 + t044;
	Real t056 = -4. * t048 * t055;
	Real t057 = t054 + t056;
	Real t058 = t044 * v2x;
	Real t059 = t046 + t058;
	Real t060 = Square(v2x);
	Real t061 = Square(x0);
	Real t062 = t051 * t061 + t011 * t048 * y0;
	Real t064 = Square(t048);
	Real t065 = t055 + t061;
	Real t066 = t048 + t060;
	Real t067 = 2. * t066;
	Real t068 = 1. + t055 + (-2. + x0) * x0;
	Real t069 = t048 + t050;
	Real t070 = 2. * t069;	
	Real t072 = 1. / t043;
	Real t073 = t043 + t061;
	Real t076 = Sqrt(t069);
	Real t077 = 1. / t076;
	Real t083 = Square(v1x);
	Real t084 = Square(v1y);
	Real t085 = t083 + t084;
	Real t086 = Square(yy0);
	Real t087 = Square(yz0);
	Real t088 = t087 + Square(yx0);
	Real t089 = t086 + t088;
	Real t090 = t018 + t019 + t085 + t089;
	Real t091 = Sqrt(t090);
	Real t092 = 1. + t015 + t089;
	Real t093 = Sqrt(t092);
	Real t094 = Sqrt(t089);
	Real t095 = Square(t014);
	Real t096 = t083 + t095;
	Real t098 = Sqrt(t096);
	Real t097 = 1. / (t096 * t098);
	Real t100 = Sqrt(t085);
	Real t099 = 1. / (t085 * t100);
	Real t114 = Sqrt(t066);
	Real t106 = 1. / (t066 * t114);
	Real t107 = 1. / (t069 * t076);
	Real t110 = t006 * t050;
	Real t141 = Square(t005);
	Real t111 = t005 * t141;
	Real t113 = 1. / t048;
	Real t115 = 1. / t114;
	Real t144 = Square(t008);
	Real t116 = t008 * t144;
	Real t117 = v2x * t060;
	Real t118 = Sqrt(t043 * t048);
	Real t119 = 1. / t118;
	Real t120 = t008 + t065 + t066;
	Real t121 = Sqrt(t120);
	Real t122 = t005 + t068 + t069;
	Real t123 = Sqrt(t122);
	Real t124 = Sqrt(t065);
	Real t125 = Sqrt(t068);
	Real t127 = t048 * t073;
	Real t128 = t003 * t011 + t055 * t060 + t127;
	Real t130 = -(t002 * t009) + t050 * t055 + t048 * (1. + t010 + t073);
	Real t132 = 1. / t066;
	Real t133 = 1. / t069;
	Real t134 = t086 + t087;
	Real t135 = t014 * t016 * t018 + t088 * t095 + t083 * (1. + t015 + t134);
	Real t137 = t017 * t019 + t084 * t088 + t083 * t134;
	Real t139 = 1. / t085;
	Real t140 = 1. / t096;
	Real t153 = -t003 + t066;
	Real t154 = t153 - t024 * v2x + x0;
	Real t156 = -t001 + t060 + t025 * v2y;
	Real t168 = 2. * t114 * t124;
	Real t169 = t043 * t051 * y0;
	Real t170 = t043 * v2y;
	Real t171 = t048 * t065;
	Real t172 = 2. * t171;
	Real t173 = -(t066 * x0);
	Real t174 = -(t067 * x0);
	Real t175 = 4. * t117;
	Real t176 = -(t175 * x0);
	Real t177 = t033 * t066 + t176;
	Real t180 = -4. * t061 * t066;
	Real t182 = t116 * v2x;
	Real t183 = ArcTan(t003 * t119);
	Real t184 = 4. * t118 * t183;
	Real t185 = t031 - 2. * t173;
	Real t186 = t076 * t123;
	Real t187 = t060 * v2y;
	Real t188 = t043 * t064;
	Real t190 = t043 * x0;
	Real t191 = t048 * t068;
	Real t192 = -2. * t191;
	Real t193 = -(t000 * t069);
	Real t194 = t036 * t069;
	Real t195 = -(t000 * t070);
	Real t196 = 4. * t052 * t069;
	Real t197 = t006 * t111;
	Real t198 = 4. * t110;
	Real t199 = -(t000 * t198);
	Real t205 = t076 * t125;
	Real t206 = t114 * t121;
	Real t207 = t003 - t066;
	Real t226 = -ArcTan(t025 * t119 * v2y);
	Real t227 = -4. * t113 * t118 * t226;
	Real t228 = t034 - 2. * t193;
	Real t229 = t009 * t113;
	Real t234 = t020 * t140;
	Real t235 = -t017 + t083;
	Real t236 = t235 + t023 * v1y;
	Real t238 = t014 * t023 + t235;

  // negative powers
	Complex t240 = Sqrt(Complex(t049));
	Complex t112 = t049 * t240;
	Complex t241 = Sqrt(Complex(t057));

	Complex t189 = t072 * t112 * x0;
	Complex t242 = t240 * y0;
	Complex t243 = -2. * t050 * t242;
	Complex t246 = -2. * t060 * t242;
	Complex t302 = -6. * t240;
	Complex t303 = -2. * t066 * t240;
	Complex t306 = -t001 + t241;
	Complex t307 = -t182 + t144 * (t173 + (-5. * t001 + t241) * v2x) - t172 * (t174 + (t011 + t241) * v2x) - t180 * (t173 + t306 * v2x) + t008 * (4. * t171 * v2x + (t177 - t060 * t302 - t303) * x0);
	Complex t308 = 1. / t240;
	Complex t247 = Sqrt(v2y * (t062 + t246 + t030 * (t059 + t241 * x0)));
	Complex t310 = conj(t247);
	Complex t312 = t127 * t240;
	Complex t313 = t169 + t312;
	Complex t314 = 1. / t310;
	Complex t315 = -(t308 * t314) * 0.5;
	Complex t317 = -(t121 * t240);
	Complex t318 = t024 * t242;
	Complex t320 = t189 + t048 * t061 * t240;
	Complex t321 = t320 + t003 * t240 * x0;
	Complex t322 = t002 + t241;
	Complex t323 = t070 * t240;
	Complex t324 = -t197 + t141 * (t193 + t006 * (5. * t002 + t241)) + t192 * (t195 + t006 * (t035 + t241)) + t196 * (t193 + t006 * t322) + t005 * (4. * t006 * t191 + t000 * (t194 + t199 - t050 * t302 + t323));
	Complex t328 = -(t125 * t240);
	Complex t329 = t052 * t240;
	Complex t330 = t002 * t003 * t240;
	Complex t331 = t043 + t240;

	Complex t244 = Sqrt(v2y * (t053 + t243 + t006 * (t047 + t000 * t241) * v2y));
	Complex t338 = conj(t244);
	Complex t339 = 1. / t338;
	Complex t340 = t187 * (-t170 + t242);

  // logs
  Complex z1n = (-4. * (-t188 + t313 + t310 * t317 + t340 + t032 * (t170 - t240 * (v2y + y0))));
  Complex z1d = (t060 + t207 + t240);
	Complex z1 = Log(z1n/z1d);
  Complex z2 = lognegconj(z1);

  Complex z3n = -4. * (t169 + t048 * (t043 * (-t002 + t240) + t329) + t330 + t328 * t338);
  Complex z3d = (t003 + t240);
	Complex z3 = Log(z3n/z3d);
  Complex z4 = lognegconj(z3);

  Complex z5n = (-4. * (-t112 + t169 - t188 - t048 * t190 + t321 + t317 * t338 + t340 + t030 * (-t318 + v2y * (t190 + t331 - t240 * x0))));
  Complex z5d = (t003 + t050 - t069 + t240);
	Complex z5 = Log(z5n/z5d);
  Complex z6 = lognegconj(z5);

  Complex z7n = -4. * (t032 * (t170 - t242) - t124 * t240 * t310 + t313);
  Complex z7d = (t003 + t240);
	Complex z7 = Log(z7n/z7d);
  Complex z8 = lognegconj(z7);

  Complex z1c = t307 * t315;
  Complex z2c = conj(z1c);
  Complex z3c = 0.5 * t308 * t324 * t339;
  Complex z4c = conj(z3c);
  Complex z5c = z3c;
  Complex z6c = z4c;
  Complex z7c = z1c;
  Complex z8c = z2c;

	Real t362 = Log((-t017 + t085 + t091 * t098 + yy0 - v1y * (1. + yy0)) / (t016 - t017 - t021 + t093 * t098 + v1y));
	Real t363 = Log((t022 - t085 - t091 * t100) / (t022 - t094 * t100));
	Real t367 = Log((t076 * t121 + t154) / (t004 + t205));
	Real t368 = Log((t001 - t206 + t207) / (t007 - t114 * t124));
	Real t370 = Log(t008 + 2. * (t066 + t206));
	Real t371 = Log(t005 + 2. * t205);
	Real t372 = Log(t121 + v2x - x0);
	Real t374 = Log(t005 + 2. * (t069 + t186));
	Real t375 = Log(-t000 + t006 + t123);
	Real t376 = Log(t008 + t168);
	Real t380 = 2. * t115 * t376;
  
  Real t369 = Log(t008 - t011 + t048 + t055);
	Real t373 = Log(t005 - t035 + t048 + t055);
	Real t377 = Log(t055);

  Real t383_1 = t184 + (-t008 + t011) * t377 + t185 * t380;
	Real t383 = t383_1 + real(z8c * z8) + real(z7c * z7);
	Real t384 = 4. * t375;
	Real t385_1 = (-4. + t227 +  t229 * t373 - 2. * t077 * t113 * t228 * t374 + t384) * t048;
  Real t385 = t385_1 + real(z5c * z5) + real(z6c * z6);

  // S1 correction
  Real S1corr;
  if((imag(z1c) >= 0.) == (imag(z5c) >= 0.))
    S1corr = ((imag(z1) + imag(z2) - imag(z7) - imag(z8)) + (imag(z5) + imag(z6) - imag(z3) - imag(z4)));
  else
    S1corr = ((imag(z1) + imag(z2) - imag(z7) - imag(z8)) - (imag(z5) + imag(z6) - imag(z3) - imag(z4)));
  S1corr *= imag(z1c);

  * S0 = 0.5;
  * S1 = t113 * S1corr * 0.25 + (4. + 4. * t113 * t118 * t226 + t113 * (real(z1c * z1) + real(z2c * z2)) - t229 * t369 + 2. * t113 * t115 * t185 * t370 - 4. * t372) * 0.25 - (t113 * (-4. * t118 * t183 + real(z3c * z3) + real(z4c * z4) + 2. * t077 * (-t034 + t195) * t371 + t009 * t377)) * 0.25 - (t113 * t383) * 0.25 + (t113 * t385) * 0.25;
  * G0x = (-1. - v2x + 3. * x0) * ONESIXTH;
  * G0y = (-1. - v1y + 3. * yy0) * ONESIXTH;
  * G1x = (t007 * t124 * t132 + t004 * t125 * t133 - t121 * t133 * t154 + t121 * t132 * t156 - t107 * t130 * t367 + t106 * t128 * t368) * 0.5;
  * G1y = (t022 * t094 * t139 - t093 * t234 + t091 * t139 * t236 - t091 * t140 * t238 - t097 * t135 * t362 + t099 * t137 * t363) * 0.5;
}

void G3(Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0)
{
	Real u000 = Square(y0) + Square(z0);
	Real u002 = -1. + v2x;
	Real u003 = v2y * y0;
	Real u004 = v2x * x0;
	Real u005 = u002 - u003 - u004 + x0;
	Real u007 = Square(v2y);
	Real u009 = u003 + u004;
	Real u011 = u007 + Square(v2x);
	Real u012 = -2. * u003;
	Real u014 = Square(v1x);
	Real u015 = Sqrt(u014 + Square(-1. + v1y));
	Real u016 = Square(yz0);
	Real u017 = u016 + Square(yy0);
	Real u018 = u017 + Square(yx0);
	Real u019 = Sqrt(1. + u018 - 2. * yy0);
	Real u020 = u014 + Square(v1y);
	Real u021 = v1y * yy0;
	Real u022 = v1x * yx0;
	Real u023 = Sqrt(u018 + u020 - 2. * u021 - 2. * u022);
	Real u024 = -u022 + yy0;
	Real u025 = Log((u020 + u015 * u023 + u024 - v1y * (1. + yy0)) / (-1. + u015 * u019 - u021 + u024 + v1y)) / u015;
	Real u026 = Sqrt(u020);
	Real u027 = Sqrt(u018);
	Real u028 = -u021 - u022;
	Real u029 = -(Log(-((u020 + u023 * u026 + u028) / (u021 + u022 - u026 * u027))) / u026);
	Real u030 = Sqrt(u011);
	Real u031 = u000 + Square(x0);
	Real u032 = Sqrt(u031);
	Real u033 = Sqrt(-2. * u004 + u011 + u012 + u031);
	Real u034 = -(Log((u009 - u011 - u030 * u033) / (u009 - u030 * u032)) / u030);
	Real u035 = Sqrt(1. + u011 - 2. * v2x);
	Real u036 = Sqrt(1. + u031 - 2. * x0);
	Real u037 = Log((-u003 + u011 + u033 * u035 + x0 - v2x * (1. + x0)) / (u005 + u035 * u036)) / u035;

  * G3x = u034 + u037;
  * G3y = u025 + u029;
}

void S3_G3(Real *S3, Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0)
{
	Real u000 = Square(y0) + Square(z0);
	Real u001 = 1. + u000 + (-2. + x0) * x0;
	Real u002 = -1. + v2x;
	Real u003 = v2y * y0;
	Real u004 = v2x * x0;
	Real u005 = u002 - u003 - u004 + x0;
	Real u006 = 2. * u005;
	Real u007 = Square(v2y);
	Real u008 = Square(u002) + u007;
	Real u009 = u003 + u004;
	Real u010 = -2. * u009;
	Real u011 = u007 + Square(v2x);
	Real u012 = -2. * u003;
	Real u013 = 1. - x0;
	Real u014 = Square(v1x);
	Real u015 = Sqrt(u014 + Square(-1. + v1y));
	Real u016 = Square(yz0);
	Real u017 = u016 + Square(yy0);
	Real u018 = u017 + Square(yx0);
	Real u019 = Sqrt(1. + u018 - 2. * yy0);
	Real u020 = u014 + Square(v1y);
	Real u021 = v1y * yy0;
	Real u022 = v1x * yx0;
	Real u023 = Sqrt(u018 + u020 - 2. * u021 - 2. * u022);
	Real u024 = -u022 + yy0;
	Real u025 = Log((u020 + u015 * u023 + u024 - v1y * (1. + yy0)) / (-1. + u015 * u019 - u021 + u024 + v1y)) / u015;
	Real u026 = Sqrt(u020);
	Real u027 = Sqrt(u018);
	Real u028 = -u021 - u022;
	Real u029 = -(Log(-((u020 + u023 * u026 + u028) / (u021 + u022 - u026 * u027))) / u026);
	Real u030 = Sqrt(u011);
	Real u031 = u000 + Square(x0);
	Real u032 = Sqrt(u031);
	Real u033 = Sqrt(-2. * u004 + u011 + u012 + u031);
	Real u034 = -(Log((u009 - u011 - u030 * u033) / (u009 - u030 * u032)) / u030);
	Real u035 = Sqrt(1. + u011 - 2. * v2x);
	Real u036 = Sqrt(1. + u031 - 2. * x0);
	Real u037 = Log((-u003 + u011 + u033 * u035 + x0 - v2x * (1. + x0)) / (u005 + u035 * u036)) / u035;
	Real u038 = Square(u012);
	Real u039 = Sqrt(2.);
	Real u040 = -4. * u000 * u007;
	Real u041 = u038 + u040;
	Real u044 = 2. * u007;
	Real u050 = 2. * u011;
	Real u054 = u010 + u050;
	Real u056 = u038 * u054;
	Real u057 = u032 * u044;
	Real u063 = Sqrt(u001 + u006 + u008);
	Real u066 = 2. * u008;
	Real u068 = u006 + u066;
	Real u070 = u038 * u068;
	Real u071 = u036 * u044;
	Complex u042 = Sqrt(Complex(u041));
	Complex u043 = u012 - u042;
	Complex u045 = -(u012 * u042);
  Complex u069 = u045 * u068;
	Complex u046 = u038 + u045;
	Complex u047 = -2. * u000 * u007 + u046;
	Complex u048 = u010 * u042;
	Complex u049 = Sqrt(u011 * u047 + u007 * (-(u010 * u012) + u031 * u044 + u048));
	Complex u051 = (u043 + u044) / u043;
	Complex u052 = u040 + u046;
	Complex u053 = 2. * u031 * u042;
	Complex u055 = u045 * u054;
	Complex u058 = u039 * u049;
	Complex u059 = -u012 + u042;
	Complex u060 = 2. * u032 * u042;
	Complex u061 = 2. * u033 * u042;
  Complex u062log = Log((u051 * (u010 * u052 + (u057 + u058) * u060)) / (u044 * (u048 + u053 - 2. * u000 * u054) + u055 + u056 + u058 * u061));
  Complex u062 = -(u059 * v2x - u044 * x0) / u049 * u062log;
	Complex u064 = u006 * u042;
	Complex u065 = Sqrt(u008 * u047 + u007 * (-(u006 * u012) + u001 * u044 + u064));
	Complex u067 = 2. * u001 * u042;
	Complex u072 = u039 * u065;
	Complex u073 = 2. * u036 * u042;
	Complex u074 = 2. * u042 * u063;
  Complex u075log = Log((u051 * (u006 * u052 + (u071 + u072) * u073)) / (u044 * (u064 + u067 - 2. * u000 * u068) + u069 + u070 + u072 * u074));
	Complex u075 = (u013 * u044 + u002 * u059) / u065 * u075log;
	Complex u076 = u012 + u042;
	Complex u078 = conj(u049);
	Complex u083 = -(u076 * v2x + 2. * u007 * x0) / u078 * conj(u062log);
	Complex u084 = conj(u065);
  Complex u086 = (-2. * u007 * u013 + u002 * u076) / u084 * conj(u075log);

  * S3 = real((u062 + u075 + u083 + u086) / (u039 * u042));
  * G3x = u034 + u037;
  * G3y = u025 + u029;
}

void C3(Real *C3z, Real v2x, Real v2y, Real x0, Real y0, Real z0)
{
	Real t000 = -1. + v2x;
	Real t001 = v2y * y0;
	Real t002 = v2x * x0;
	Real t003 = -t001 - t002;
	Real t004 = t003 + x0;
	Real t005 = t000 + t004;
	Real t006 = 2. * t005;
	Real t007 = t001 + t002;
	Real t008 = -2. * t007;
	Real t009 = -2. * t001;
	Real t010 = 1. - x0;
	Real t011 = -(t010 * x0);
	Real t012 = t002 * y0;
	Real t017 = t004 * v2y;
	Real t018 = t009 * t017;
	Real t019 = -(t007 * v2y);
	Real t020 = Square(v2y);
	Real t021 = 1. / t020;
	Real t022 = Square(t009);
	Real t023 = Square(z0);
	Real t024 = Square(y0);
	Real t025 = t023 + t024;
	Real t026 = -4. * t020 * t025;
	Real t027 = t022 + t026;
	Real t030 = 2. * t020;
	Real t036 = 1. + t025 + (-2. + x0) * x0;
	Real t037 = Square(t000) + t020;
	Real t038 = Sqrt(t006 + t036 + t037);
	Real t039 = Square(x0);
	Real t040 = t025 + t039;
	Real t041 = Square(t020);
	Real t042 = 2. * t040 * t041;
	Real t043 = Square(v2x);
	Real t044 = t020 + t043;
	Real t045 = 2. * t044;
	Real t055 = 2. * t037;
	Real t064 = Sqrt(2.);
	Real t065 = Sqrt(1. + t044 - 2. * v2x);
	Real t066 = 1. / t065;

	Real t067 = Sqrt(1. + t040 - 2. * x0);
	Real t068 = Sqrt(t044);
	Real t069 = 1. / t068;

	Real t070 = Sqrt(-2. * t002 + t009 + t040 + t044);
	Real t071 = Sqrt(t040);
	Real t082 = t068 * t070;
	Real t083 = t065 * t067;
	Real t084 = 2. * t017 * t066;
	Real t086 = t030 * t067;
	Real t089 = 2. * t025;
	Real t091 = 2. * (t011 + t025) * t041 * y0;
	Real t092 = t011 * v2y + 2. * t024 * v2y + t000 * x0 * y0;
	Real t093 = t017 * t089 + t009 * t092;
	Real t101 = 2. * t019 * t069;
	Real t102 = t030 * t071;
	Real t108 = t012 + (2. * t024 + t039) * v2y;
	Real t111 = t006 + t055;
	Real t114 = t022 * t111;
	Real t118 = t008 + t045;
	Real t122 = t022 * t118;
	Real t136 = t084 * Log((t006 + 2. * (t037 + t038 * t065)) / (t006 + 2. * t083));
	Real t139 = t101 * Log((t008 + 2. * t068 * t071) / (t008 + 2. * (t044 + t082)));
	Complex t028 = Sqrt(Complex(t027));
	Complex t029 = t009 + t028;
	Complex t031 = t029 + t030;
	Complex t032 = 1. / t031;
	Complex t046 = Sqrt(t042 + t009 * t029 * t044 - t020 * (t008 * t029 + t025 * t045));
	Complex t047 = 1. / t046;
	Complex t048 = -(t009 * t028);
	Complex t053 = Sqrt(Complex(-(t020 * t023)));
	Complex t054 = 1. / Sqrt(v2y * (-2. * t012 * t020 + t039 * Power(v2y,3.) + v2x * v2y * ((-t023 + t024) * v2x - 2. * t053 * x0) + 2. * t043 * t053 * y0));
	Complex t056 = Sqrt(t009 * t029 * t037 + 2. * t036 * t041 - t020 * (t006 * t029 + t025 * t055));
	Complex t057 = 1. / t056;
	Complex t058 = 1. / t029;
	Complex t060 = conj(t056);
	Complex t061 = 1. / t060;
	Complex t063 = 1. / t028;
	Complex t088 = 2. * t028 * t067;
	Complex t090 = -t009 + t028;
	Complex t094 = t063 * t064;
	Complex t095 = t061 * (t018 * t090 + t091 + t020 * (-(t028 * t092) + t093)) * t094;
	Complex t097 = t056 * t064;
	Complex t098 = t027 + t009 * t028;
	Complex t099 = t057 * (t018 * t029 - t091 - t020 * (t028 * t092 + t093)) * t094;
	Complex t100 = t058 * (t088 * (t086 + t097) - t006 * t098);
	Complex t104 = 2. * t028 * t071;
	Complex t105 = t054 * v2y * x0 * (t001 * v2x + t053 * v2x - t020 * x0);
	Complex t107 = t046 * t064;
	Complex t109 = t047 * t094 * (t009 * t019 * t029 - t020 * (t019 * t089 + t009 * t108 + t028 * t108) - t042 * y0);
	Complex t110 = t058 * (-(t008 * t098) + t104 * (t102 + t107));
	Complex t112 = t048 * t111;
	Complex t113 = 2. * t028 * t036;
	Complex t116 = t028 + t089;
	Complex t117 = 2. * t028 * t038 * t097 + t112 - t114 + t030 * (4. * t025 * t037 + t113 + t006 * t116);
	Complex t119 = t048 * t118;
	Complex t120 = 2. * t028 * t040;
	Complex t121 = 2. * t028 * t070;
	Complex t124 = t119 + t030 * (4. * t025 * t044 + t008 * t116 + t120) + t107 * t121 - t122;

  Complex t128log = Log((t031 * t110) / t124);
  Complex t134log = Log(t032 * t117 / t100);

	Complex t134 = t099 * t134log;
	Complex t135 = t095 * conj(t134log);
	Complex t137 = t109 * t128log;
	Complex t138 = t105 * conj(t128log);

  * C3z = real(0.5 * t021 * (t134 + t135 + t136 + t137 + t138 + t139));
}

void S3_G3_C3(Real *S3, Real *G3x, Real *G3y, Real *C3z, Real v2x, Real v2y, Real x0, Real y0, Real z0, Real v1x, Real v1y, Real yx0, Real yy0, Real yz0)
{
	Real t000 = -1. + v2x;
	Real t001 = v2y * y0;
	Real t002 = v2x * x0;
	Real t003 = -t001 - t002;
	Real t004 = t003 + x0;
	Real t005 = t000 + t004;
	Real t006 = 2. * t005;
	Real t007 = t001 + t002;
	Real t008 = -2. * t007;
	Real t009 = -2. * t001;
	Real t010 = 1. - x0;
	Real t011 = -(t010 * x0);
	Real t012 = t002 * y0;
	Real t013 = v1y * yy0;
	Real t014 = v1x * yx0;
	Real t015 = -t014 + yy0;
	Real t016 = -t013 - t014;
	Real t017 = t004 * v2y;
	Real t018 = t009 * t017;
	Real t019 = -(t007 * v2y);
	Real t020 = Square(v2y);
	Real t021 = 1. / t020;
	Real t022 = Square(t009);
	Real t023 = Square(z0);
	Real t024 = Square(y0);
	Real t025 = t023 + t024;
	Real t026 = -4. * t020 * t025;
	Real t027 = t022 + t026;
	Real t030 = 2. * t020;
	Real t036 = 1. + t025 + (-2. + x0) * x0;
	Real t037 = Square(t000) + t020;
	Real t038 = Sqrt(t006 + t036 + t037);
	Real t039 = Square(x0);
	Real t040 = t025 + t039;
	Real t041 = Square(t020);
	Real t042 = 2. * t040 * t041;
	Real t043 = Square(v2x);
	Real t044 = t020 + t043;
	Real t045 = 2. * t044;
	Real t055 = 2. * t037;
	Real t064 = Sqrt(2.);
	Real t065 = Sqrt(1. + t044 - 2. * v2x);
	Real t066 = 1. / t065;
	Real t067 = Sqrt(1. + t040 - 2. * x0);
	Real t068 = Sqrt(t044);
	Real t069 = 1. / t068;
	Real t070 = Sqrt(-2. * t002 + t009 + t040 + t044);
	Real t071 = Sqrt(t040);
	Real t072 = Square(yz0);
	Real t073 = t072 + Square(yy0);
	Real t074 = t073 + Square(yx0);
	Real t075 = Sqrt(t074);
	Real t076 = Square(v1x);
	Real t077 = t076 + Square(v1y);
	Real t078 = Sqrt(t077);
	Real t079 = Sqrt(-2. * t013 - 2. * t014 + t074 + t077);
	Real t080 = Sqrt(1. + t074 - 2. * yy0);
	Real t081 = Sqrt(t076 + Square(-1. + v1y));
	Real t082 = t068 * t070;
	Real t083 = t065 * t067;
	Real t084 = 2. * t017 * t066;
	Real t086 = t030 * t067;
	Real t089 = 2. * t025;
	Real t091 = 2. * (t011 + t025) * t041 * y0;
	Real t092 = t011 * v2y + 2. * t024 * v2y + t000 * x0 * y0;
	Real t093 = t017 * t089 + t009 * t092;
	Real t101 = 2. * t019 * t069;
	Real t102 = t030 * t071;
	Real t108 = t012 + (2. * t024 + t039) * v2y;
	Real t111 = t006 + t055;
	Real t114 = t022 * t111;
	Real t118 = t008 + t045;
	Real t122 = t022 * t118;
	Real t125 = -(Log(-((t016 + t077 + t078 * t079) / (t013 + t014 - t075 * t078))) / t078);
	Real t126 = Log((t015 + t077 + t079 * t081 - v1y * (1. + yy0)) / (-1. - t013 + t015 + t080 * t081 + v1y)) / t081;
	Real t143 = t066 * Log((-t001 + t044 + t065 * t070 + x0 - v2x * (1. + x0)) / (t005 + t083));
	Real t144 = -(t069 * Log((t007 - t044 - t082) / (t007 - t068 * t071)));
	Real t136 = t084 * Log((t006 + 2. * (t037 + t038 * t065)) / (t006 + 2. * t083));
	Real t139 = t101 * Log((t008 + 2. * t068 * t071) / (t008 + 2. * (t044 + t082)));
	Complex t028 = Sqrt(Complex(t027));
	Complex t029 = t009 + t028;
	Complex t031 = t029 + t030;
	Complex t032 = 1. / t031;
	Complex t046 = Sqrt(t042 + t009 * t029 * t044 - t020 * (t008 * t029 + t025 * t045));
	Complex t047 = 1. / t046;
	Complex t048 = -(t009 * t028);
  Complex t052 = conj(t046);
	Complex t053 = Sqrt(Complex(-(t020 * t023)));
	Complex t054 = 1. / Sqrt(v2y * (-2. * t012 * t020 + t039 * Power(v2y,3.) + v2x * v2y * ((-t023 + t024) * v2x - 2. * t053 * x0) + 2. * t043 * t053 * y0));
	Complex t056 = Sqrt(t009 * t029 * t037 + 2. * t036 * t041 - t020 * (t006 * t029 + t025 * t055));
	Complex t057 = 1. / t056;
	Complex t058 = 1. / t029;
	Complex t060 = conj(t056);
	Complex t061 = 1. / t060;
	Complex t063 = 1. / t028;
	Complex t088 = 2. * t028 * t067;
	Complex t090 = -t009 + t028;
	Complex t094 = t063 * t064;
	Complex t095 = t061 * (t018 * t090 + t091 + t020 * (-(t028 * t092) + t093)) * t094;
	Complex t097 = t056 * t064;
	Complex t098 = t027 + t009 * t028;
	Complex t099 = t057 * (t018 * t029 - t091 - t020 * (t028 * t092 + t093)) * t094;
	Complex t100 = t058 * (t088 * (t086 + t097) - t006 * t098);
	Complex t104 = 2. * t028 * t071;
	Complex t105 = t054 * v2y * x0 * (t001 * v2x + t053 * v2x - t020 * x0);
	Complex t107 = t046 * t064;
	Complex t109 = t047 * t094 * (t009 * t019 * t029 - t020 * (t019 * t089 + t009 * t108 + t028 * t108) - t042 * y0);
	Complex t110 = t058 * (-(t008 * t098) + t104 * (t102 + t107));
	Complex t112 = t048 * t111;
	Complex t113 = 2. * t028 * t036;
	Complex t116 = t028 + t089;
	Complex t117 = 2. * t028 * t038 * t097 + t112 - t114 + t030 * (4. * t025 * t037 + t113 + t006 * t116);
	Complex t119 = t048 * t118;
	Complex t120 = 2. * t028 * t040;
	Complex t121 = 2. * t028 * t070;
	Complex t124 = t119 + t030 * (4. * t025 * t044 + t008 * t116 + t120) + t107 * t121 - t122;
  Complex t127log = Log((t031 * t100) / t117);
	Complex t127 = (-2. * t010 * t020 + t000 * t029) * t057 * t127log;
  Complex t128log = Log((t031 * t110) / t124);
	Complex t128 = -(t047 * (t029 * v2x + 2. * t020 * x0) * t128log);
  Complex t129 = t061 * (t010 * t030 + t000 * t090) * conj(t127log);
  Complex t130 = -(((t090 * v2x - t030 * x0) * conj(t128log)) / t052);
  Complex t134log = Log(t032 * t117 / t100);
	Complex t134 = t099 * t134log;
	Complex t135 = t095 * conj(t134log);
	Complex t137 = t109 * t128log;
	Complex t138 = t105 * conj(t128log);

  * S3 = real((t063 * (t127 + t128 + t129 + t130)) / t064);
  * G3x = t143 + t144;
  * G3y = t125 + t126;
  * C3z = real(0.5 * t021 * (t134 + t135 + t136 + t137 + t138 + t139));
}

// this is incorrect at the corners
void G3overlog_edgelimit(Real *G3x, Real *G3y, Real v2x, Real v2y, Real x0, Real y0, Real v1x, Real v1y)
{
  Real y = y0 / v2y;
  Real x = x0 - v2x * y;
  
  if(fabs(x) <= fabs(y) && fabs(x) <= fabs(1. - y - x)) {
    * G3x = 1. / Sqrt(Square(v2x) + Square(v2y));
    * G3y = 0.;
  } else if(fabs(y) <= fabs(x) && fabs(y) <= fabs(1. - y - x)) {
    * G3x = 0.;
    * G3y = 1. / Sqrt(Square(v1x) + Square(v1y));
  } else {
    * G3x = -1. / Sqrt(Square(v2x - 1.) + Square(v2y));
    * G3y = -1. / Sqrt(Square(v1y - 1.) + Square(v1x));
  }
}

void C3overlog_edgelimit(Real *C3z, Real v2x, Real v2y, Real x0, Real y0)
{ 
  Real y = y0 / v2y;
  Real x = x0 - v2x * y;

  if(fabs(x) <= fabs(y) && fabs(x) <= fabs(1. - y - x)) {
    Real t000 = Square(v2y);
    * C3z = -y0 * Sqrt(Square(v2x) + t000) / t000;   
  } else if(fabs(y) <= fabs(x) && fabs(y) <= fabs(1. - y - x)) {
    * C3z = x0 / v2y;
  } else {
    Real t000 = v2x - 1.;
    Real t001 = Square(t000);
    Real t002 = Square(v2y);
    Real t003 = t001 + t002;
    Real t004 = 1. / Sqrt(t003);
    * C3z = (t000 * v2y + t003 * y0) * t004 / t002;
  }
}
