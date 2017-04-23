#ifndef TRIANGLEQUADRATURE_H
#define TRIANGLEQUADRATURE_H

#include "Vector3.h"
#include "FunctionArgs.h"
#include "MathUtils.h"

class TriangleIntegrand 
{
 public:
  virtual void weighPoint(Real w, Vector3 &v)=0;
};

class TriangleQuadrature
{
 public:
  virtual void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2)=0;
};

class TriangleQuadratureGauss3 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		v = v0 + 0.166666666666666 * v1 + 0.166666666666667 * v2; integrand.weighPoint(0.333333333333333,v);
		v = v0 + 0.166666666666666 * v1 + 0.666666666666667 * v2; integrand.weighPoint(0.333333333333333,v);
		v = v0 + 0.666666666666666 * v1 + 0.166666666666667 * v2; integrand.weighPoint(0.333333333333333,v);
	}
};

class TriangleQuadratureGauss6 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		v = v0 + 0.8168475729804 * v1 + 0.0915762135097999 * v2; integrand.weighPoint(0.1099517436553,v);
		v = v0 + 0.0915762135097 * v1 + 0.0915762135098 * v2; integrand.weighPoint(0.1099517436553,v);
		v = v0 + 0.0915762135097 * v1 + 0.8168475729805 * v2; integrand.weighPoint(0.1099517436553,v);
		v = v0 + 0.4459484909159 * v1 + 0.445948490916 * v2; integrand.weighPoint(0.223381589678,v);
		v = v0 + 0.4459484909159 * v1 + 0.1081030181681 * v2; integrand.weighPoint(0.223381589678,v);
		v = v0 + 0.108103018168 * v1 + 0.445948490916 * v2; integrand.weighPoint(0.223381589678,v);
	}
};

class TriangleQuadratureGauss10 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		//WARNING: Evaluating on boundary
		//WARNING: Evaluating on boundary
		v = v0 + v2; integrand.weighPoint(0.0131356049752,v);
		//WARNING: Evaluating on boundary
		v = v0; integrand.weighPoint(0.0131358306034,v);
		//WARNING: Evaluating on boundary
		v = v0 + v1; integrand.weighPoint(0.01370819738,v);
		v = v0 + 0.0598527250105 * v1 + 0.672819921871 * v2; integrand.weighPoint(0.11741919329115,v);
		v = v0 + 0.0598535871057 * v1 + 0.2673288599482 * v2; integrand.weighPoint(0.1174206119134,v);
		v = v0 + 0.2634233538452 * v1 + 0.6716530111494 * v2; integrand.weighPoint(0.1240125896557,v);
		v = v0 + 0.2634249770929 * v1 + 0.0649251690029 * v2; integrand.weighPoint(0.12401524612605,v);
		v = v0 + 0.6652178176747 * v1 + 0.2693789366453 * v2; integrand.weighPoint(0.12593023027645,v);
		v = v0 + 0.6652178055941 * v1 + 0.0654054874919 * v2; integrand.weighPoint(0.1259330266829,v);
		v = v0 + 0.3226461603077 * v1 + 0.3386799893027 * v2; integrand.weighPoint(0.2252894690957,v);
	}
};

class TriangleQuadratureGauss15 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		//WARNING: Evaluating on boundary
		//WARNING: Evaluating on boundary
		v = v0; integrand.weighPoint(0.0051279087046,v);
		//WARNING: Evaluating on boundary
		v = v0 + v1; integrand.weighPoint(0.0051279087046,v);
		//WARNING: Evaluating on boundary
		v = v0 + v2; integrand.weighPoint(0.0051279087046,v);
		v = v0 + 0.1738960507346 * v1 + 0.0421382841642 * v2; integrand.weighPoint(0.05580235233235,v);
		v = v0 + 0.0421382841643 * v1 + 0.7839656651012 * v2; integrand.weighPoint(0.05580235233235,v);
		v = v0 + 0.7839656651013 * v1 + 0.0421382841642 * v2; integrand.weighPoint(0.05580235233235,v);
		v = v0 + 0.7839656651013 * v1 + 0.1738960507345 * v2; integrand.weighPoint(0.05580235233235,v);
		v = v0 + 0.0421382841643 * v1 + 0.1738960507345 * v2; integrand.weighPoint(0.05580235233235,v);
		v = v0 + 0.1738960507346 * v1 + 0.7839656651012 * v2; integrand.weighPoint(0.05580235233235,v);
		v = v0 + 0.0512238276496 * v1 + 0.4743880861752 * v2; integrand.weighPoint(0.08398877976675,v);
		v = v0 + 0.4743880861751 * v1 + 0.0512238276497 * v2; integrand.weighPoint(0.08398877976675,v);
		v = v0 + 0.4743880861751 * v1 + 0.4743880861752 * v2; integrand.weighPoint(0.08398877976675,v);
		v = v0 + 0.238561530018 * v1 + 0.5228769399639 * v2; integrand.weighPoint(0.1326119401973,v);
		v = v0 + 0.238561530018 * v1 + 0.2385615300181 * v2; integrand.weighPoint(0.1326119401973,v);
		v = v0 + 0.5228769399638 * v1 + 0.2385615300181 * v2; integrand.weighPoint(0.1326119401973,v);
	}
};

class TriangleQuadratureGauss21 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		v = v0 + 0.9096219804312 * v1 + 0.0451890097844 * v2; integrand.weighPoint(0.0259935710323,v);
		v = v0 + 0.0451890097844 * v1 + 0.9096219804312 * v2; integrand.weighPoint(0.0259935710323,v);
		v = v0 + 0.0451890097844 * v1 + 0.0451890097844 * v2; integrand.weighPoint(0.0259935710323,v);
		v = v0 + 0.2220631655373 * v1 + 0.0304243617288 * v2; integrand.weighPoint(0.0353517050892,v);
		v = v0 + 0.7475124727339 * v1 + 0.0304243617288 * v2; integrand.weighPoint(0.0353517050892,v);
		v = v0 + 0.0304243617288 * v1 + 0.2220631655373 * v2; integrand.weighPoint(0.0353517050892,v);
		v = v0 + 0.0304243617288 * v1 + 0.7475124727339 * v2; integrand.weighPoint(0.0353517050892,v);
		v = v0 + 0.2220631655373 * v1 + 0.7475124727339 * v2; integrand.weighPoint(0.0353517050892,v);
		v = v0 + 0.7475124727339 * v1 + 0.2220631655373 * v2; integrand.weighPoint(0.0353517050892,v);
		v = v0 + 0.6447187277637 * v1 + 0.2182900709714 * v2; integrand.weighPoint(0.0454695380476,v);
		v = v0 + 0.1369912012649 * v1 + 0.2182900709714 * v2; integrand.weighPoint(0.0454695380476,v);
		v = v0 + 0.2182900709714 * v1 + 0.6447187277637 * v2; integrand.weighPoint(0.0454695380476,v);
		v = v0 + 0.1369912012649 * v1 + 0.6447187277637 * v2; integrand.weighPoint(0.0454695380476,v);
		v = v0 + 0.6447187277637 * v1 + 0.1369912012649 * v2; integrand.weighPoint(0.0454695380476,v);
		v = v0 + 0.2182900709714 * v1 + 0.1369912012649 * v2; integrand.weighPoint(0.0454695380476,v);
		v = v0 + 0.4815198347833 * v1 + 0.4815198347833 * v2; integrand.weighPoint(0.051617202569,v);
		v = v0 + 0.4815198347833 * v1 + 0.0369603304334 * v2; integrand.weighPoint(0.051617202569,v);
		v = v0 + 0.0369603304333999 * v1 + 0.4815198347833 * v2; integrand.weighPoint(0.051617202569,v);
		v = v0 + 0.403603979818 * v1 + 0.1927920403641 * v2; integrand.weighPoint(0.09408007345835,v);
		v = v0 + 0.1927920403642 * v1 + 0.4036039798179 * v2; integrand.weighPoint(0.09408007345835,v);
		v = v0 + 0.403603979818 * v1 + 0.4036039798179 * v2; integrand.weighPoint(0.09408007345835,v);
	}
};

class TriangleQuadratureGauss28 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		//WARNING: Evaluating on boundary
		//WARNING: Evaluating on boundary
		v = v0 + 0.0548295549826 * v1 + 0.9451704450174 * v2; integrand.weighPoint(0.00570412470165,v);
		//WARNING: Evaluating on boundary
		v = v0 + 0.0548295549827 * v1; integrand.weighPoint(0.00570412470165,v);
		v = v0 + 0.00254917970569998 * v1 + 0.0685505797224 * v2; integrand.weighPoint(0.006634564286,v);
		v = v0 + 0.00254917970590007 * v1 + 0.9289002405717 * v2; integrand.weighPoint(0.006634564286,v);
		v = v0 + 0.9513463288769 * v1 + 0.0243268355616 * v2; integrand.weighPoint(0.0077932886675,v);
		v = v0 + 0.8442498415177 * v1 + 0.0277838749488 * v2; integrand.weighPoint(0.0204137390214,v);
		v = v0 + 0.8442498415175 * v1 + 0.1279662835337 * v2; integrand.weighPoint(0.02041373902145,v);
		v = v0 + 0.2214568982983 * v1 + 0.7498347588657 * v2; integrand.weighPoint(0.0289924832558,v);
		v = v0 + 0.2214568982984 * v1 + 0.028708342836 * v2; integrand.weighPoint(0.0289924832558,v);
		v = v0 + 0.0274390027907999 * v1 + 0.2497602062385 * v2; integrand.weighPoint(0.03006926238315,v);
		v = v0 + 0.0274390027907 * v1 + 0.7228007909707 * v2; integrand.weighPoint(0.03006926238315,v);
		v = v0 + 0.0808923150164 * v1 + 0.8325513856997 * v2; integrand.weighPoint(0.03126369442165,v);
		v = v0 + 0.0808923150163 * v1 + 0.0865562992839 * v2; integrand.weighPoint(0.03126369442165,v);
		v = v0 + 0.6634854224837 * v1 + 0.0303526617491 * v2; integrand.weighPoint(0.0319842160752,v);
		v = v0 + 0.6634854224834 * v1 + 0.3061619157675 * v2; integrand.weighPoint(0.0319842160752,v);
		v = v0 + 0.0262778809906 * v1 + 0.4868610595047 * v2; integrand.weighPoint(0.03306629360805,v);
		v = v0 + 0.1576639552764 * v1 + 0.1765456154219 * v2; integrand.weighPoint(0.033425161841,v);
		v = v0 + 0.1576639552763 * v1 + 0.6657904293016 * v2; integrand.weighPoint(0.03342516184105,v);
		v = v0 + 0.4411221503971 * v1 + 0.5295657488669 * v2; integrand.weighPoint(0.03434521529885,v);
		v = v0 + 0.4411221503973 * v1 + 0.029312100736 * v2; integrand.weighPoint(0.03434521529885,v);
		v = v0 + 0.7110652351218 * v1 + 0.1444673824391 * v2; integrand.weighPoint(0.05013587719295,v);
		v = v0 + 0.1338444159539 * v1 + 0.536181572905 * v2; integrand.weighPoint(0.05715683920495,v);
		v = v0 + 0.1338444159539 * v1 + 0.3299740111409 * v2; integrand.weighPoint(0.05715683920495,v);
		v = v0 + 0.3050701621215 * v1 + 0.1437790861923 * v2; integrand.weighPoint(0.0611824073376,v);
		v = v0 + 0.3050701621215 * v1 + 0.5511507516862 * v2; integrand.weighPoint(0.0611824073376,v);
		v = v0 + 0.5122313975512 * v1 + 0.1529619437161 * v2; integrand.weighPoint(0.0697211167089,v);
		v = v0 + 0.5122313975512 * v1 + 0.3348066587327 * v2; integrand.weighPoint(0.0697211167089,v);
		v = v0 + 0.3139633003706 * v1 + 0.3430183498147 * v2; integrand.weighPoint(0.0872188914591,v);
	}
};

class TriangleQuadratureGauss36 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		v = v0 + 0.0264005354564 * v1 + 0.9493059293846 * v2; integrand.weighPoint(0.00831204993785,v);
		v = v0 + 0.9492111441638 * v1 + 0.024269513064 * v2; integrand.weighPoint(0.0083405849889,v);
		v = v0 + 0.0242806010012 * v1 + 0.0265067966437 * v2; integrand.weighPoint(0.00834152845335,v);
		v = v0 + 0.5198907823888 * v1 + 0.4767316412363 * v2; integrand.weighPoint(0.00878404350415,v);
		v = v0 + 0.00434058727969999 * v1 + 0.5198921829102 * v2; integrand.weighPoint(0.00922373309225,v);
		v = v0 + 0.4753304100327 * v1 + 0.0055912706202 * v2; integrand.weighPoint(0.0098971205094,v);
		v = v0 + 0.1249164206061 * v1 + 0.0133996048618 * v2; integrand.weighPoint(0.01017701979275,v);
		v = v0 + 0.013773591874 * v1 + 0.8613054321334 * v2; integrand.weighPoint(0.010342643197,v);
		v = v0 + 0.8613700828781 * v1 + 0.1247733717358 * v2; integrand.weighPoint(0.0104135683043,v);
		v = v0 + 0.1349674584555 * v1 + 0.8438438351223 * v2; integrand.weighPoint(0.01589098891395,v);
		v = v0 + 0.0213139566951 * v1 + 0.135456364583 * v2; integrand.weighPoint(0.01602360176205,v);
		v = v0 + 0.8432285381479 * v1 + 0.0213482820656 * v2; integrand.weighPoint(0.0160303840573,v);
		v = v0 + 0.6689226826307 * v1 + 0.0221919663014 * v2; integrand.weighPoint(0.02153829795915,v);
		v = v0 + 0.0225929525442 * v1 + 0.3089012879389 * v2; integrand.weighPoint(0.02192367076695,v);
		v = v0 + 0.3081745044122 * v1 + 0.6691709943321 * v2; integrand.weighPoint(0.02196048363665,v);
		v = v0 + 0.0266766436122 * v1 + 0.6924718155106 * v2; integrand.weighPoint(0.02399759618455,v);
		v = v0 + 0.2808829905923 * v1 + 0.0268723345026 * v2; integrand.weighPoint(0.02419031303665,v);
		v = v0 + 0.6921288579659 * v1 + 0.2810093973222 * v2; integrand.weighPoint(0.02424337116875,v);
		v = v0 + 0.0884640100944001 * v1 + 0.7973581413586 * v2; integrand.weighPoint(0.0278482244012,v);
		v = v0 + 0.1145385569148 * v1 + 0.0879806508791 * v2; integrand.weighPoint(0.0280513182178,v);
		v = v0 + 0.7962172144978 * v1 + 0.1145020561128 * v2; integrand.weighPoint(0.02825950618465,v);
		v = v0 + 0.2260607987623 * v1 + 0.6686904119922 * v2; integrand.weighPoint(0.0344644945335,v);
		v = v0 + 0.1061926087428 * v1 + 0.2275051631832 * v2; integrand.weighPoint(0.03586066680445,v);
		v = v0 + 0.6637623701232 * v1 + 0.1054572561221 * v2; integrand.weighPoint(0.0363726960488,v);
		v = v0 + 0.3120876443802 * v1 + 0.5174064398658 * v2; integrand.weighPoint(0.03944036683685,v);
		v = v0 + 0.1742882171748 * v1 + 0.3170523855209 * v2; integrand.weighPoint(0.0405057172756,v);
		v = v0 + 0.504746977606 * v1 + 0.1810706361659 * v2; integrand.weighPoint(0.04128626495275,v);
		v = v0 + 0.0703944642332 * v1 + 0.4678594539804 * v2; integrand.weighPoint(0.0421022283665,v);
		v = v0 + 0.4684056461834 * v1 + 0.4622856042085 * v2; integrand.weighPoint(0.04217927666525,v);
		v = v0 + 0.4623686935063 * v1 + 0.0724357805669 * v2; integrand.weighPoint(0.0425984934244,v);
		v = v0 + 0.128997910293 * v1 + 0.6131395039177 * v2; integrand.weighPoint(0.0451422664026,v);
		v = v0 + 0.2587011398612 * v1 + 0.1300360834609 * v2; integrand.weighPoint(0.04571415717425,v);
		v = v0 + 0.6113104035182 * v1 + 0.2581713828884 * v2; integrand.weighPoint(0.04581395327045,v);
		v = v0 + 0.3356556038355 * v1 + 0.2362005969817 * v2; integrand.weighPoint(0.0512786687448,v);
		v = v0 + 0.2331977907682 * v1 + 0.4311026308588 * v2; integrand.weighPoint(0.05165798307065,v);
		v = v0 + 0.4238561751788 * v1 + 0.3456013949376 * v2; integrand.weighPoint(0.05179271835965,v);
	}
};

class TriangleQuadratureGauss45 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		//WARNING: Evaluating on boundary
		//WARNING: Evaluating on boundary
		v = v0 + v2; integrand.weighPoint(0.0005308355995,v);
		//WARNING: Evaluating on boundary
		v = v0; integrand.weighPoint(0.0005308355995,v);
		//WARNING: Evaluating on boundary
		v = v0 + v1; integrand.weighPoint(0.0005308355995,v);
		v = v0 + 0.927528685716 * v1 + 0.0151382269814 * v2; integrand.weighPoint(0.00657301180505,v);
		v = v0 + 0.0151382269814 * v1 + 0.927528685716 * v2; integrand.weighPoint(0.00657301180505,v);
		v = v0 + 0.0151382269814 * v1 + 0.0573330873026 * v2; integrand.weighPoint(0.00657301180505,v);
		v = v0 + 0.927528685716 * v1 + 0.0573330873026 * v2; integrand.weighPoint(0.00657301180505,v);
		v = v0 + 0.0573330873026 * v1 + 0.0151382269814 * v2; integrand.weighPoint(0.00657301180505,v);
		v = v0 + 0.0573330873026 * v1 + 0.927528685716 * v2; integrand.weighPoint(0.00657301180505,v);
		v = v0 + 0.0180654989724 * v1 + 0.1659719969565 * v2; integrand.weighPoint(0.01214409634745,v);
		v = v0 + 0.1659719969565 * v1 + 0.0180654989724 * v2; integrand.weighPoint(0.01214409634745,v);
		v = v0 + 0.0180654989724 * v1 + 0.8159625040711 * v2; integrand.weighPoint(0.01214409634745,v);
		v = v0 + 0.1659719969565 * v1 + 0.8159625040711 * v2; integrand.weighPoint(0.01214409634745,v);
		v = v0 + 0.8159625040711 * v1 + 0.0180654989724 * v2; integrand.weighPoint(0.01214409634745,v);
		v = v0 + 0.8159625040711 * v1 + 0.1659719969565 * v2; integrand.weighPoint(0.01214409634745,v);
		v = v0 + 0.6647637544849 * v1 + 0.0186886898773 * v2; integrand.weighPoint(0.0158399933166,v);
		v = v0 + 0.3165475556378 * v1 + 0.0186886898773 * v2; integrand.weighPoint(0.0158399933166,v);
		v = v0 + 0.3165475556378 * v1 + 0.6647637544849 * v2; integrand.weighPoint(0.0158399933166,v);
		v = v0 + 0.6647637544849 * v1 + 0.3165475556378 * v2; integrand.weighPoint(0.0158399933166,v);
		v = v0 + 0.0186886898773 * v1 + 0.6647637544849 * v2; integrand.weighPoint(0.0158399933166,v);
		v = v0 + 0.0186886898773 * v1 + 0.3165475556378 * v2; integrand.weighPoint(0.0158399933166,v);
		v = v0 + 0.4903668903754 * v1 + 0.4903668903754 * v2; integrand.weighPoint(0.0174658973518,v);
		v = v0 + 0.4903668903754 * v1 + 0.0192662192492 * v2; integrand.weighPoint(0.0174658973518,v);
		v = v0 + 0.0192662192492001 * v1 + 0.4903668903754 * v2; integrand.weighPoint(0.0174658973518,v);
		v = v0 + 0.0875134669581999 * v1 + 0.8249730660837 * v2; integrand.weighPoint(0.01918322669725,v);
		v = v0 + 0.8249730660838 * v1 + 0.0875134669581 * v2; integrand.weighPoint(0.01918322669725,v);
		v = v0 + 0.0875134669582 * v1 + 0.0875134669581 * v2; integrand.weighPoint(0.01918322669725,v);
		v = v0 + 0.6984608540614 * v1 + 0.2079865423167 * v2; integrand.weighPoint(0.0289184745605,v);
		v = v0 + 0.2079865423168 * v1 + 0.6984608540613 * v2; integrand.weighPoint(0.0289184745605,v);
		v = v0 + 0.6984608540614 * v1 + 0.0935526036219 * v2; integrand.weighPoint(0.0289184745605,v);
		v = v0 + 0.2079865423168 * v1 + 0.0935526036219 * v2; integrand.weighPoint(0.0289184745605,v);
		v = v0 + 0.093552603622 * v1 + 0.2079865423167 * v2; integrand.weighPoint(0.0289184745605,v);
		v = v0 + 0.093552603622 * v1 + 0.6984608540613 * v2; integrand.weighPoint(0.0289184745605,v);
		v = v0 + 0.3645018421384 * v1 + 0.5380088595149 * v2; integrand.weighPoint(0.0362910843697,v);
		v = v0 + 0.538008859515 * v1 + 0.0974892983467 * v2; integrand.weighPoint(0.0362910843697,v);
		v = v0 + 0.3645018421384 * v1 + 0.0974892983467 * v2; integrand.weighPoint(0.0362910843697,v);
		v = v0 + 0.0974892983468 * v1 + 0.3645018421383 * v2; integrand.weighPoint(0.0362910843697,v);
		v = v0 + 0.0974892983468001 * v1 + 0.5380088595149 * v2; integrand.weighPoint(0.0362910843697,v);
		v = v0 + 0.538008859515 * v1 + 0.3645018421383 * v2; integrand.weighPoint(0.0362910843697,v);
		v = v0 + 0.2217145894874 * v1 + 0.5565708210253 * v2; integrand.weighPoint(0.04489282620535,v);
		v = v0 + 0.2217145894874 * v1 + 0.2217145894873 * v2; integrand.weighPoint(0.04489282620535,v);
		v = v0 + 0.5565708210254 * v1 + 0.2217145894873 * v2; integrand.weighPoint(0.04489282620535,v);
		v = v0 + 0.3860471669296 * v1 + 0.2279056661408 * v2; integrand.weighPoint(0.05172722668085,v);
		v = v0 + 0.3860471669296 * v1 + 0.3860471669296 * v2; integrand.weighPoint(0.05172722668085,v);
		v = v0 + 0.2279056661408 * v1 + 0.3860471669296 * v2; integrand.weighPoint(0.05172722668085,v);
	}
};

class TriangleQuadratureGauss55 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		//WARNING: Evaluating on boundary
		//WARNING: Evaluating on boundary
		v = v0; integrand.weighPoint(0.00031012999255,v);
		//WARNING: Evaluating on boundary
		v = v0 + v2; integrand.weighPoint(0.0003157587356,v);
		//WARNING: Evaluating on boundary
		v = v0 + v1; integrand.weighPoint(0.00035433007795,v);
		v = v0 + 0.0551287671789 * v1 + 0.0049848744634 * v2; integrand.weighPoint(0.0027581858084,v);
		v = v0 + 0.00697876983250001 * v1 + 0.9386405618617 * v2; integrand.weighPoint(0.0031346203828,v);
		v = v0 + 0.9379635488139 * v1 + 0.0526424462697 * v2; integrand.weighPoint(0.0039265704413,v);
		v = v0 + 0.0366619396287 * v1 + 0.9469035517351 * v2; integrand.weighPoint(0.0047275741932,v);
		v = v0 + 0.0167139052971 * v1 + 0.0363373677167 * v2; integrand.weighPoint(0.00489122556355,v);
		v = v0 + 0.9422171452433 * v1 + 0.0151224541799 * v2; integrand.weighPoint(0.00499308217445,v);
		v = v0 + 0.1183956993897 * v1 + 0.8693773510664 * v2; integrand.weighPoint(0.0068776909408,v);
		v = v0 + 0.0121386193179 * v1 + 0.1204917285774 * v2; integrand.weighPoint(0.007048958902,v);
		v = v0 + 0.1385492010741 * v1 + 0.015776396787 * v2; integrand.weighPoint(0.00748234321685,v);
		v = v0 + 0.0156119497522 * v1 + 0.8448120870375 * v2; integrand.weighPoint(0.0078048751806,v);
		v = v0 + 0.8547168651185 * v1 + 0.0135009605584 * v2; integrand.weighPoint(0.0078841846674,v);
		v = v0 + 0.8386769935164 * v1 + 0.1455274938536 * v2; integrand.weighPoint(0.00878972731915,v);
		v = v0 + 0.2478839574656 * v1 + 0.0155697540908 * v2; integrand.weighPoint(0.0102056920135,v);
		v = v0 + 0.248047467522 * v1 + 0.737983689445 * v2; integrand.weighPoint(0.0104781439308,v);
		v = v0 + 0.015448912419 * v1 + 0.7297615689771 * v2; integrand.weighPoint(0.0105356706499,v);
		v = v0 + 0.0140536794130001 * v1 + 0.2543076683315 * v2; integrand.weighPoint(0.0108823380101,v);
		v = v0 + 0.7146506475258 * v1 + 0.2696239795791 * v2; integrand.weighPoint(0.01111442043495,v);
		v = v0 + 0.7192913200045 * v1 + 0.0144783956308 * v2; integrand.weighPoint(0.0112093346841,v);
		v = v0 + 0.0734816524385999 * v1 + 0.05916794104 * v2; integrand.weighPoint(0.01150613084965,v);
		v = v0 + 0.0623723757982 * v1 + 0.8634782575061 * v2; integrand.weighPoint(0.011840695125,v);
		v = v0 + 0.5649475096402 * v1 + 0.4191238955238 * v2; integrand.weighPoint(0.0128732321684,v);
		v = v0 + 0.4034716050786 * v1 + 0.5809222921146 * v2; integrand.weighPoint(0.0128978400804,v);
		v = v0 + 0.3930653729865 * v1 + 0.0159251452651 * v2; integrand.weighPoint(0.0129036163805,v);
		v = v0 + 0.0158528135007 * v1 + 0.5806700368104 * v2; integrand.weighPoint(0.01301716160295,v);
		v = v0 + 0.0155759225172 * v1 + 0.4149495146302 * v2; integrand.weighPoint(0.01328840708045,v);
		v = v0 + 0.8560287620759 * v1 + 0.0761218678591 * v2; integrand.weighPoint(0.01328923809155,v);
		v = v0 + 0.5576521717416 * v1 + 0.0157509692312 * v2; integrand.weighPoint(0.0133766164619,v);
		v = v0 + 0.1587119179689 * v1 + 0.7741898312421 * v2; integrand.weighPoint(0.01878939033205,v);
		v = v0 + 0.1652570272881 * v1 + 0.0819119495639 * v2; integrand.weighPoint(0.01915329470975,v);
		v = v0 + 0.0669143759151 * v1 + 0.1577128457292 * v2; integrand.weighPoint(0.01924248475125,v);
		v = v0 + 0.0806983742471 * v1 + 0.7503943099742 * v2; integrand.weighPoint(0.0194809912926,v);
		v = v0 + 0.7604352659813 * v1 + 0.0708311507268 * v2; integrand.weighPoint(0.01973020557735,v);
		v = v0 + 0.7415758664793 * v1 + 0.1762996626771 * v2; integrand.weighPoint(0.0206182389049,v);
		v = v0 + 0.2903549683338 * v1 + 0.0807744953317 * v2; integrand.weighPoint(0.02564362192415,v);
		v = v0 + 0.6134213394958 * v1 + 0.3054373589776 * v2; integrand.weighPoint(0.02582028209675,v);
		v = v0 + 0.0803401946049001 * v1 + 0.6227485988871 * v2; integrand.weighPoint(0.02591150211345,v);
		v = v0 + 0.2985210536283 * v1 + 0.6247247149546 * v2; integrand.weighPoint(0.02642639940905,v);
		v = v0 + 0.0765491844989 * v1 + 0.3011485821166 * v2; integrand.weighPoint(0.02692527865135,v);
		v = v0 + 0.611711534687 * v1 + 0.0779098365079 * v2; integrand.weighPoint(0.02709476646595,v);
		v = v0 + 0.4577148746462 * v1 + 0.4603633038351 * v2; integrand.weighPoint(0.0292368573222,v);
		v = v0 + 0.446142332819 * v1 + 0.0821554006797 * v2; integrand.weighPoint(0.02964315841815,v);
		v = v0 + 0.081583155086 * v1 + 0.463756503389 * v2; integrand.weighPoint(0.02971791383745,v);
		v = v0 + 0.1876630852575 * v1 + 0.6422277808188 * v2; integrand.weighPoint(0.03159001279315,v);
		v = v0 + 0.1695702133257 * v1 + 0.1898293537256 * v2; integrand.weighPoint(0.03164634225765,v);
		v = v0 + 0.634777673094 * v1 + 0.1739955685343 * v2; integrand.weighPoint(0.0320353680886,v);
		v = v0 + 0.3315770162524 * v1 + 0.4798914070406 * v2; integrand.weighPoint(0.0406020297959,v);
		v = v0 + 0.187871344419 * v1 + 0.3348356598119 * v2; integrand.weighPoint(0.0407218756765,v);
		v = v0 + 0.1915053180981 * v1 + 0.4957972197259 * v2; integrand.weighPoint(0.04073396006205,v);
		v = v0 + 0.311122038515 * v1 + 0.1927553668904 * v2; integrand.weighPoint(0.0407525274042,v);
		v = v0 + 0.4910178879872 * v1 + 0.3161015807261 * v2; integrand.weighPoint(0.04075823324695,v);
		v = v0 + 0.4745065744894 * v1 + 0.189489280129 * v2; integrand.weighPoint(0.04084655298115,v);
		v = v0 + 0.3319148427341 * v1 + 0.3343571021811 * v2; integrand.weighPoint(0.04616091672655,v);
	}
};

class TriangleQuadratureGauss91 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		v = v0 + 0.9928951216156 * v1 + 0.00355243919220001 * v2; integrand.weighPoint(0.00033522182195,v);
		v = v0 + 0.00355243919219994 * v1 + 0.9928951216156 * v2; integrand.weighPoint(0.00033522182195,v);
		v = v0 + 0.00355243919219999 * v1 + 0.0035524391922 * v2; integrand.weighPoint(0.00033522182195,v);
		v = v0 + 0.0358552797177 * v1 + 0.0087898929093 * v2; integrand.weighPoint(0.0022736304037,v);
		v = v0 + 0.955354827373 * v1 + 0.00878989290929999 * v2; integrand.weighPoint(0.0022736304037,v);
		v = v0 + 0.00878989290929996 * v1 + 0.0358552797177 * v2; integrand.weighPoint(0.0022736304037,v);
		v = v0 + 0.955354827373 * v1 + 0.0358552797177 * v2; integrand.weighPoint(0.0022736304037,v);
		v = v0 + 0.0358552797177 * v1 + 0.955354827373 * v2; integrand.weighPoint(0.0022736304037,v);
		v = v0 + 0.00878989290929999 * v1 + 0.955354827373 * v2; integrand.weighPoint(0.0022736304037,v);
		v = v0 + 0.00524053759359999 * v1 + 0.1082329745017 * v2; integrand.weighPoint(0.002603879266,v);
		v = v0 + 0.1082329745018 * v1 + 0.0052405375935 * v2; integrand.weighPoint(0.002603879266,v);
		v = v0 + 0.8865264879048 * v1 + 0.1082329745017 * v2; integrand.weighPoint(0.002603879266,v);
		v = v0 + 0.1082329745018 * v1 + 0.8865264879047 * v2; integrand.weighPoint(0.002603879266,v);
		v = v0 + 0.00524053759359999 * v1 + 0.8865264879047 * v2; integrand.weighPoint(0.002603879266,v);
		v = v0 + 0.8865264879048 * v1 + 0.0052405375935 * v2; integrand.weighPoint(0.002603879266,v);
		v = v0 + 0.046639743215 * v1 + 0.90672051357 * v2; integrand.weighPoint(0.00327177164435,v);
		v = v0 + 0.90672051357 * v1 + 0.046639743215 * v2; integrand.weighPoint(0.00327177164435,v);
		v = v0 + 0.0466397432149999 * v1 + 0.046639743215 * v2; integrand.weighPoint(0.00327177164435,v);
		v = v0 + 0.784152030177 * v1 + 0.00827592412839998 * v2; integrand.weighPoint(0.00463689207665,v);
		v = v0 + 0.00827592412840006 * v1 + 0.784152030177 * v2; integrand.weighPoint(0.00463689207665,v);
		v = v0 + 0.00827592412840003 * v1 + 0.2075720456946 * v2; integrand.weighPoint(0.00463689207665,v);
		v = v0 + 0.2075720456946 * v1 + 0.784152030177 * v2; integrand.weighPoint(0.00463689207665,v);
		v = v0 + 0.784152030177 * v1 + 0.2075720456946 * v2; integrand.weighPoint(0.00463689207665,v);
		v = v0 + 0.2075720456946 * v1 + 0.0082759241284 * v2; integrand.weighPoint(0.00463689207665,v);
		v = v0 + 0.8827043562574 * v1 + 0.0314836947701 * v2; integrand.weighPoint(0.00479688913115,v);
		v = v0 + 0.0858119489725001 * v1 + 0.0314836947701 * v2; integrand.weighPoint(0.00479688913115,v);
		v = v0 + 0.8827043562574 * v1 + 0.0858119489725 * v2; integrand.weighPoint(0.00479688913115,v);
		v = v0 + 0.0314836947701 * v1 + 0.8827043562574 * v2; integrand.weighPoint(0.00479688913115,v);
		v = v0 + 0.0314836947701 * v1 + 0.0858119489725 * v2; integrand.weighPoint(0.00479688913115,v);
		v = v0 + 0.0858119489725 * v1 + 0.8827043562574 * v2; integrand.weighPoint(0.00479688913115,v);
		v = v0 + 0.3216071005549 * v1 + 0.0095150760625 * v2; integrand.weighPoint(0.00571239045835,v);
		v = v0 + 0.6688778233825 * v1 + 0.321607100555 * v2; integrand.weighPoint(0.00571239045835,v);
		v = v0 + 0.3216071005549 * v1 + 0.6688778233826 * v2; integrand.weighPoint(0.00571239045835,v);
		v = v0 + 0.00951507606239993 * v1 + 0.321607100555 * v2; integrand.weighPoint(0.00571239045835,v);
		v = v0 + 0.00951507606239987 * v1 + 0.6688778233826 * v2; integrand.weighPoint(0.00571239045835,v);
		v = v0 + 0.6688778233825 * v1 + 0.00951507606250002 * v2; integrand.weighPoint(0.00571239045835,v);
		v = v0 + 0.5520140671206 * v1 + 0.00998597856810002 * v2; integrand.weighPoint(0.0058608482087,v);
		v = v0 + 0.4379999543113 * v1 + 0.5520140671206 * v2; integrand.weighPoint(0.0058608482087,v);
		v = v0 + 0.00998597856810002 * v1 + 0.5520140671206 * v2; integrand.weighPoint(0.0058608482087,v);
		v = v0 + 0.5520140671206 * v1 + 0.4379999543113 * v2; integrand.weighPoint(0.0058608482087,v);
		v = v0 + 0.00998597856810002 * v1 + 0.4379999543113 * v2; integrand.weighPoint(0.0058608482087,v);
		v = v0 + 0.4379999543113 * v1 + 0.0099859785681 * v2; integrand.weighPoint(0.0058608482087,v);
		v = v0 + 0.1619974933733 * v1 + 0.0405093994119 * v2; integrand.weighPoint(0.0094098577616,v);
		v = v0 + 0.7974931072147 * v1 + 0.1619974933734 * v2; integrand.weighPoint(0.0094098577616,v);
		v = v0 + 0.1619974933733 * v1 + 0.7974931072148 * v2; integrand.weighPoint(0.0094098577616,v);
		v = v0 + 0.0405093994118 * v1 + 0.7974931072148 * v2; integrand.weighPoint(0.0094098577616,v);
		v = v0 + 0.0405093994118 * v1 + 0.1619974933734 * v2; integrand.weighPoint(0.0094098577616,v);
		v = v0 + 0.7974931072147 * v1 + 0.0405093994119 * v2; integrand.weighPoint(0.0094098577616,v);
		v = v0 + 0.227156889609 * v1 + 0.3864215551955 * v2; integrand.weighPoint(0.01176304901355,v);
		v = v0 + 0.3864215551955 * v1 + 0.227156889609 * v2; integrand.weighPoint(0.01176304901355,v);
		v = v0 + 0.3864215551955 * v1 + 0.3864215551955 * v2; integrand.weighPoint(0.01176304901355,v);
		v = v0 + 0.0954935310335 * v1 + 0.0954935310336 * v2; integrand.weighPoint(0.01177857330755,v);
		v = v0 + 0.0954935310335 * v1 + 0.8090129379329 * v2; integrand.weighPoint(0.01177857330755,v);
		v = v0 + 0.8090129379328 * v1 + 0.0954935310336 * v2; integrand.weighPoint(0.01177857330755,v);
		v = v0 + 0.6774734280561 * v1 + 0.0479840480721 * v2; integrand.weighPoint(0.0134123103715,v);
		v = v0 + 0.2745425238718 * v1 + 0.6774734280561 * v2; integrand.weighPoint(0.0134123103715,v);
		v = v0 + 0.2745425238718 * v1 + 0.0479840480721 * v2; integrand.weighPoint(0.0134123103715,v);
		v = v0 + 0.0479840480721 * v1 + 0.2745425238718 * v2; integrand.weighPoint(0.0134123103715,v);
		v = v0 + 0.0479840480721 * v1 + 0.6774734280561 * v2; integrand.weighPoint(0.0134123103715,v);
		v = v0 + 0.6774734280561 * v1 + 0.2745425238718 * v2; integrand.weighPoint(0.0134123103715,v);
		v = v0 + 0.0516677930989 * v1 + 0.5429849622344 * v2; integrand.weighPoint(0.01571448883895,v);
		v = v0 + 0.5429849622344 * v1 + 0.4053472446667 * v2; integrand.weighPoint(0.01571448883895,v);
		v = v0 + 0.5429849622344 * v1 + 0.0516677930989 * v2; integrand.weighPoint(0.01571448883895,v);
		v = v0 + 0.4053472446667 * v1 + 0.0516677930989 * v2; integrand.weighPoint(0.01571448883895,v);
		v = v0 + 0.4053472446667 * v1 + 0.5429849622344 * v2; integrand.weighPoint(0.01571448883895,v);
		v = v0 + 0.0516677930989 * v1 + 0.4053472446667 * v2; integrand.weighPoint(0.01571448883895,v);
		v = v0 + 0.7054113116873 * v1 + 0.1068148267588 * v2; integrand.weighPoint(0.01685980960795,v);
		v = v0 + 0.1068148267589 * v1 + 0.1877738615539 * v2; integrand.weighPoint(0.01685980960795,v);
		v = v0 + 0.187773861554 * v1 + 0.1068148267588 * v2; integrand.weighPoint(0.01685980960795,v);
		v = v0 + 0.187773861554 * v1 + 0.7054113116872 * v2; integrand.weighPoint(0.01685980960795,v);
		v = v0 + 0.1068148267589 * v1 + 0.7054113116872 * v2; integrand.weighPoint(0.01685980960795,v);
		v = v0 + 0.7054113116873 * v1 + 0.1877738615539 * v2; integrand.weighPoint(0.01685980960795,v);
		v = v0 + 0.5747817297348 * v1 + 0.3057122990643 * v2; integrand.weighPoint(0.02138726471065,v);
		v = v0 + 0.3057122990643 * v1 + 0.5747817297348 * v2; integrand.weighPoint(0.02138726471065,v);
		v = v0 + 0.3057122990643 * v1 + 0.1195059712009 * v2; integrand.weighPoint(0.02138726471065,v);
		v = v0 + 0.1195059712009 * v1 + 0.3057122990643 * v2; integrand.weighPoint(0.02138726471065,v);
		v = v0 + 0.1195059712009 * v1 + 0.5747817297348 * v2; integrand.weighPoint(0.02138726471065,v);
		v = v0 + 0.5747817297348 * v1 + 0.1195059712009 * v2; integrand.weighPoint(0.02138726471065,v);
		v = v0 + 0.2009377128318 * v1 + 0.2009377128319 * v2; integrand.weighPoint(0.02205694663685,v);
		v = v0 + 0.2009377128318 * v1 + 0.5981245743363 * v2; integrand.weighPoint(0.02205694663685,v);
		v = v0 + 0.5981245743362 * v1 + 0.2009377128319 * v2; integrand.weighPoint(0.02205694663685,v);
		v = v0 + 0.4717864543322 * v1 + 0.3121360256673 * v2; integrand.weighPoint(0.0230734797342,v);
		v = v0 + 0.4717864543322 * v1 + 0.2160775200005 * v2; integrand.weighPoint(0.0230734797342,v);
		v = v0 + 0.3121360256674 * v1 + 0.4717864543321 * v2; integrand.weighPoint(0.0230734797342,v);
		v = v0 + 0.2160775200006 * v1 + 0.4717864543321 * v2; integrand.weighPoint(0.0230734797342,v);
		v = v0 + 0.2160775200006 * v1 + 0.3121360256673 * v2; integrand.weighPoint(0.0230734797342,v);
		v = v0 + 0.3121360256674 * v1 + 0.2160775200005 * v2; integrand.weighPoint(0.0230734797342,v);
		v = v0 + 0.1246840192302 * v1 + 0.4376579903849 * v2; integrand.weighPoint(0.0234576234312,v);
		v = v0 + 0.4376579903848 * v1 + 0.1246840192303 * v2; integrand.weighPoint(0.0234576234312,v);
		v = v0 + 0.4376579903848 * v1 + 0.4376579903849 * v2; integrand.weighPoint(0.0234576234312,v);
		v = v0 + 0.3333333333334 * v1 + 0.3333333333333 * v2; integrand.weighPoint(0.02755999901735,v);
	}
};

class TriangleQuadratureGauss105 : public TriangleQuadrature
{
 public:
	void integrate(TriangleIntegrand &integrand, Vector3 &v0, Vector3 &v1, Vector3 &v2) {
		Vector3 v = v0;
		v = v0 + 0.000851425939200068 * v1 + 0.9903676436772 * v2; integrand.weighPoint(0.00032191491305,v);
		v = v0 + 0.000851546954800041 * v1 + 0.0087809216232 * v2; integrand.weighPoint(0.0003219206538,v);
		v = v0 + 0.9637056319111 * v1 + 0.0335914404439 * v2; integrand.weighPoint(0.0005067367855,v);
		v = v0 + 0.9637061838766 * v1 + 0.00270289467100001 * v2; integrand.weighPoint(0.0005067376288,v);
		v = v0 + 0.9816648578343 * v1 + 0.00916763530509997 * v2; integrand.weighPoint(0.00098399649675,v);
		v = v0 + 0.0239694640786 * v1 + 0.0084737176656 * v2; integrand.weighPoint(0.0016733656892,v);
		v = v0 + 0.0239693363967001 * v1 + 0.9675569435345 * v2; integrand.weighPoint(0.0016733669604,v);
		v = v0 + 0.9244433107346 * v1 + 0.0676784943862 * v2; integrand.weighPoint(0.00214366616875,v);
		v = v0 + 0.9244432863009 * v1 + 0.00787816592909999 * v2; integrand.weighPoint(0.00214367299425,v);
		v = v0 + 0.00867585037660002 * v1 + 0.0442974541187 * v2; integrand.weighPoint(0.0021501900686,v);
		v = v0 + 0.00867585678329996 * v1 + 0.9470266676487 * v2; integrand.weighPoint(0.0021501924549,v);
		v = v0 + 0.0774021329986 * v1 + 0.0081735455132 * v2; integrand.weighPoint(0.00284673146025,v);
		v = v0 + 0.077402034151 * v1 + 0.9144244234031 * v2; integrand.weighPoint(0.0028467320067,v);
		v = v0 + 0.3669315272539 * v1 + 0.383323243472 * v2; integrand.weighPoint(0.00308219340075,v);
		v = v0 + 0.366931608594 * v1 + 0.2497451268005 * v2; integrand.weighPoint(0.0030822378209,v);
		v = v0 + 0.00878208369970002 * v1 + 0.1035328809446 * v2; integrand.weighPoint(0.00310072567955,v);
		v = v0 + 0.00878208398630009 * v1 + 0.887684993184 * v2; integrand.weighPoint(0.0031007265976,v);
		v = v0 + 0.8519553084408 * v1 + 0.1403190991974 * v2; integrand.weighPoint(0.0034818165147,v);
		v = v0 + 0.8519551640269 * v1 + 0.00772559346239998 * v2; integrand.weighPoint(0.0034818165921,v);
		v = v0 + 0.00857664664219998 * v1 + 0.1809642523926 * v2; integrand.weighPoint(0.003753312886,v);
		v = v0 + 0.0085766480949 * v1 + 0.8104590515334 * v2; integrand.weighPoint(0.00375331322825,v);
		v = v0 + 0.1586221111639 * v1 + 0.0083010939677 * v2; integrand.weighPoint(0.00395373841695,v);
		v = v0 + 0.1586220547482 * v1 + 0.8330768545392 * v2; integrand.weighPoint(0.00395373862425,v);
		v = v0 + 0.9303185324371 * v1 + 0.0348406969482 * v2; integrand.weighPoint(0.00401766723115,v);
		v = v0 + 0.0085730472444 * v1 + 0.7173981847948 * v2; integrand.weighPoint(0.0043981720537,v);
		v = v0 + 0.00857304708359996 * v1 + 0.2740287304386 * v2; integrand.weighPoint(0.0043981724056,v);
		v = v0 + 0.7523163959504 * v1 + 0.0081859182262 * v2; integrand.weighPoint(0.0045652097858,v);
		v = v0 + 0.7523165247478 * v1 + 0.2394975566677 * v2; integrand.weighPoint(0.00456521068055,v);
		v = v0 + 0.5087422955238 * v1 + 0.4843740892687 * v2; integrand.weighPoint(0.00464108743755,v);
		v = v0 + 0.5087422281352 * v1 + 0.0068836232949 * v2; integrand.weighPoint(0.0046410907831,v);
		v = v0 + 0.00784646977520009 * v1 + 0.4960767529507 * v2; integrand.weighPoint(0.0047249903089,v);
		v = v0 + 0.00827395325159996 * v1 + 0.3804323691239 * v2; integrand.weighPoint(0.0047313734242,v);
		v = v0 + 0.00827395531219999 * v1 + 0.6112936466533 * v2; integrand.weighPoint(0.0047313742647,v);
		v = v0 + 0.2612122106775 * v1 + 0.0083987179701 * v2; integrand.weighPoint(0.00477778861425,v);
		v = v0 + 0.2612121935954 * v1 + 0.7303890895407 * v2; integrand.weighPoint(0.00477778964215,v);
		v = v0 + 0.3795998344693 * v1 + 0.0075475979695 * v2; integrand.weighPoint(0.0048069421244,v);
		v = v0 + 0.3795998554381 * v1 + 0.6128525484582 * v2; integrand.weighPoint(0.0048069423413,v);
		v = v0 + 0.6360700856766 * v1 + 0.3559773826721 * v2; integrand.weighPoint(0.0049995762106,v);
		v = v0 + 0.6360699771038 * v1 + 0.00795253585019999 * v2; integrand.weighPoint(0.0049995775925,v);
		v = v0 + 0.0452529356689 * v1 + 0.0437233665345 * v2; integrand.weighPoint(0.00501506596385,v);
		v = v0 + 0.0452529587388 * v1 + 0.9110236807446 * v2; integrand.weighPoint(0.0050150673318,v);
		v = v0 + 0.8644489029883 * v1 + 0.0967030908282 * v2; integrand.weighPoint(0.00624683380925,v);
		v = v0 + 0.8644487939678 * v1 + 0.0388479942386 * v2; integrand.weighPoint(0.00624683630625,v);
		v = v0 + 0.8253546468297 * v1 + 0.0873226620391 * v2; integrand.weighPoint(0.00700986545685,v);
		v = v0 + 0.1092937008808 * v1 + 0.8485617789108 * v2; integrand.weighPoint(0.0071668108448,v);
		v = v0 + 0.1092936604124 * v1 + 0.0421445420915 * v2; integrand.weighPoint(0.00716681360625,v);
		v = v0 + 0.0454642723664 * v1 + 0.1067435942472 * v2; integrand.weighPoint(0.007680207137,v);
		v = v0 + 0.0454642782456001 * v1 + 0.8477921328146 * v2; integrand.weighPoint(0.00768020917125,v);
		v = v0 + 0.7749692956401 * v1 + 0.0416340521608 * v2; integrand.weighPoint(0.0092261912807,v);
		v = v0 + 0.7749694261903 * v1 + 0.183396519693 * v2; integrand.weighPoint(0.0092261931573,v);
		v = v0 + 0.0446768545588 * v1 + 0.1941599202852 * v2; integrand.weighPoint(0.00979169917865,v);
		v = v0 + 0.0446768591918 * v1 + 0.7611632153938 * v2; integrand.weighPoint(0.0097917009997,v);
		v = v0 + 0.1980794644241 * v1 + 0.0439826608586 * v2; integrand.weighPoint(0.0098816375671,v);
		v = v0 + 0.1980795245297 * v1 + 0.7579378242308 * v2; integrand.weighPoint(0.00988163833385,v);
		v = v0 + 0.4267053387646 * v1 + 0.5363186076436 * v2; integrand.weighPoint(0.00994031955095,v);
		v = v0 + 0.4267052084723 * v1 + 0.0369760780935 * v2; integrand.weighPoint(0.0099403242888,v);
		v = v0 + 0.1086475957534 * v1 + 0.7912267093545 * v2; integrand.weighPoint(0.0103590919242,v);
		v = v0 + 0.1086475751803 * v1 + 0.1001257554673 * v2; integrand.weighPoint(0.01035909674465,v);
		v = v0 + 0.5462720157265 * v1 + 0.4157413128558 * v2; integrand.weighPoint(0.010447153572,v);
		v = v0 + 0.54627189095 * v1 + 0.0379867061535 * v2; integrand.weighPoint(0.0104471625978,v);
		v = v0 + 0.3072752281824 * v1 + 0.0420141226713 * v2; integrand.weighPoint(0.01074322869425,v);
		v = v0 + 0.3072753221478 * v1 + 0.6507105645084 * v2; integrand.weighPoint(0.01074322930035,v);
		v = v0 + 0.6653825532262 * v1 + 0.2920626023484 * v2; integrand.weighPoint(0.0111109066518,v);
		v = v0 + 0.6653824346007 * v1 + 0.0425548546753 * v2; integrand.weighPoint(0.01111090801015,v);
		v = v0 + 0.0417238992814999 * v1 + 0.4193031469005 * v2; integrand.weighPoint(0.01116726527275,v);
		v = v0 + 0.0417239077901 * v1 + 0.538972909361 * v2; integrand.weighPoint(0.01116726893695,v);
		v = v0 + 0.0443175354138 * v1 + 0.3007352636162 * v2; integrand.weighPoint(0.0112379462473,v);
		v = v0 + 0.0443175396352 * v1 + 0.6549471812731 * v2; integrand.weighPoint(0.011237949022,v);
		v = v0 + 0.2793619097663 * v1 + 0.3453980130752 * v2; integrand.weighPoint(0.01148506979225,v);
		v = v0 + 0.2793619021541 * v1 + 0.3752400695673 * v2; integrand.weighPoint(0.0114851697219,v);
		v = v0 + 0.7407159136052 * v1 + 0.1598308695187 * v2; integrand.weighPoint(0.0116399188051,v);
		v = v0 + 0.7407158680283 * v1 + 0.0994531960132 * v2; integrand.weighPoint(0.0116399213753,v);
		v = v0 + 0.1078087907409 * v1 + 0.7124585430924 * v2; integrand.weighPoint(0.01347415998235,v);
		v = v0 + 0.1078087815817 * v1 + 0.179732772224 * v2; integrand.weighPoint(0.01347416535535,v);
		v = v0 + 0.1932232537189 * v1 + 0.7001701784175 * v2; integrand.weighPoint(0.0140219379005,v);
		v = v0 + 0.1932232240227 * v1 + 0.1066065855677 * v2; integrand.weighPoint(0.01402193823035,v);
		v = v0 + 0.2941048385403 * v1 + 0.6065647984796 * v2; integrand.weighPoint(0.0143763135086,v);
		v = v0 + 0.294104805071 * v1 + 0.0993303896769 * v2; integrand.weighPoint(0.01437631936355,v);
		v = v0 + 0.6443394877768 * v1 + 0.2533381579528 * v2; integrand.weighPoint(0.01494904145315,v);
		v = v0 + 0.6443393848873 * v1 + 0.1023223826189 * v2; integrand.weighPoint(0.01494904613795,v);
		v = v0 + 0.1064271224208 * v1 + 0.2769502060575 * v2; integrand.weighPoint(0.0154502179258,v);
		v = v0 + 0.1064271406267 * v1 + 0.6166227900624 * v2; integrand.weighPoint(0.0154502192978,v);
		v = v0 + 0.4114292791126 * v1 + 0.4981522637001 * v2; integrand.weighPoint(0.0157015508544,v);
		v = v0 + 0.4114292187603 * v1 + 0.0904185045149 * v2; integrand.weighPoint(0.01570155369775,v);
		v = v0 + 0.5333349622924 * v1 + 0.3738418516908 * v2; integrand.weighPoint(0.0159595776512,v);
		v = v0 + 0.5333348715981 * v1 + 0.092823258479 * v2; integrand.weighPoint(0.0159595834189,v);
		v = v0 + 0.4956640233896 * v1 + 0.2521680925697 * v2; integrand.weighPoint(0.0160714962031,v);
		v = v0 + 0.1006919236962 * v1 + 0.390558054433 * v2; integrand.weighPoint(0.0165197800694,v);
		v = v0 + 0.1006919445608 * v1 + 0.5087501437661 * v2; integrand.weighPoint(0.01651978159145,v);
		v = v0 + 0.302712049135 * v1 + 0.5266738039554 * v2; integrand.weighPoint(0.01780845477945,v);
		v = v0 + 0.3027119981151 * v1 + 0.1706142257537 * v2; integrand.weighPoint(0.0178084638027,v);
		v = v0 + 0.3924363387485 * v1 + 0.2588055084886 * v2; integrand.weighPoint(0.0182870594999,v);
		v = v0 + 0.392436291228 * v1 + 0.3487583491703 * v2; integrand.weighPoint(0.0182870757602,v);
		v = v0 + 0.5289863257983 * v1 + 0.3013522183964 * v2; integrand.weighPoint(0.0182988823495,v);
		v = v0 + 0.5289862229906 * v1 + 0.1696615963219 * v2; integrand.weighPoint(0.01829890269445,v);
		v = v0 + 0.2835055815763 * v1 + 0.4584741774478 * v2; integrand.weighPoint(0.0184972840057,v);
		v = v0 + 0.2835055320791 * v1 + 0.2580203819011 * v2; integrand.weighPoint(0.01849728875295,v);
		v = v0 + 0.6302202611951 * v1 + 0.1848898704551 * v2; integrand.weighPoint(0.01870268118935,v);
		v = v0 + 0.1947647667466 * v1 + 0.1921611994069 * v2; integrand.weighPoint(0.01877751291585,v);
		v = v0 + 0.1947647850617 * v1 + 0.6130740398389 * v2; integrand.weighPoint(0.0187775156265,v);
		v = v0 + 0.4168845502985 * v1 + 0.1650613336416 * v2; integrand.weighPoint(0.0194443846743,v);
		v = v0 + 0.416884615872 * v1 + 0.4180541199244 * v2; integrand.weighPoint(0.0194443854171,v);
		v = v0 + 0.1858075255146 * v1 + 0.2982719005229 * v2; integrand.weighPoint(0.0196352821774,v);
		v = v0 + 0.1858075529888 * v1 + 0.5159205534362 * v2; integrand.weighPoint(0.01963529012585,v);
		v = v0 + 0.1802211079868 * v1 + 0.4098894317792 * v2; integrand.weighPoint(0.01993834399155,v);
	}
};

class Integrand
{
 public:
  virtual void weighPoint(Real w, Vector3 &v)=0;
};

class Quadrature
{
 public:
  virtual void integrate(Integrand &integrand, Vector3 &v0, Vector3 &v1)=0;
};

class QuadratureGauss4 : public Quadrature
{
 public:
  void integrate(Integrand &integrand, Vector3 &v0, Vector3 &v1) {
    Vector3 v = v0 + 0.06943184420297371238802675 * v1;
    integrand.weighPoint(0.17392742256872692868653195,v);
    v = v0 + 0.3300094782075718675986671 * v1;
    integrand.weighPoint(0.32607257743127307131346805,v);
    v = v0 + 0.6699905217924281324013329 * v1;
    integrand.weighPoint(0.32607257743127307131346805,v);
    v = v0 + 0.93056815579702628761197325 * v1;
    integrand.weighPoint(0.17392742256872692868653195,v);
  }
};

class QuadratureGauss7 : public Quadrature
{
 public:
  void integrate(Integrand &integrand, Vector3 &v0, Vector3 &v1) {
    Vector3 v = v0 + 0.02544604382862073773690515 * v1;
    integrand.weighPoint(0.0647424830844348466353057,v);
    v = v0 + 0.1292344072003027800680676 * v1;
    integrand.weighPoint(0.1398526957446383339507339,v);
    v = v0 + 0.2970774243113014165466968 * v1;
    integrand.weighPoint(0.1909150252525594724751849,v);
    v = v0 + 0.5 * v1;
    integrand.weighPoint(0.208979591836734693877551,v);
    v = v0 + 0.7029225756886985834533032 * v1;
    integrand.weighPoint(0.1909150252525594724751849,v);
    v = v0 + 0.8707655927996972199319324 * v1;
    integrand.weighPoint(0.1398526957446383339507339,v);
    v = v0 + 0.97455395617137926226309485 * v1;
    integrand.weighPoint(0.0647424830844348466353057,v);
  }
};


class QuadratureGauss10 : public Quadrature
{
 public:
  void integrate(Integrand &integrand, Vector3 &v0, Vector3 &v1) {
    Vector3 v = v0 + 0.013046735741414139961018 * v1;
    integrand.weighPoint(0.0333356721543440687967844,v);
    v = v0 + 0.06746831665550774463395165 * v1;
    integrand.weighPoint(0.07472567457529029657288815,v);
    v = v0 + 0.1602952158504877968828363 * v1;
    integrand.weighPoint(0.10954318125799102199776745,v);
    v = v0 + 0.28330230293537640460036705 * v1;
    integrand.weighPoint(0.13463335965499817754561345,v);
    v = v0 + 0.425562830509184394557587 * v1;
    integrand.weighPoint(0.1477621123573764350869465,v);
    v = v0 + 0.574437169490815605442413 * v1;
    integrand.weighPoint(0.1477621123573764350869465,v);
    v = v0 + 0.71669769706462359539963295 * v1;
    integrand.weighPoint(0.13463335965499817754561345,v);
    v = v0 + 0.8397047841495122031171637 * v1;
    integrand.weighPoint(0.10954318125799102199776745,v);
    v = v0 + 0.93253168334449225536604835 * v1;
    integrand.weighPoint(0.07472567457529029657288815,v);
    v = v0 + 0.986953264258585860038982 * v1;
    integrand.weighPoint(0.0333356721543440687967844,v);
  }
};



class QuadratureGauss15 : public Quadrature
{
 public:
  void integrate(Integrand &integrand, Vector3 &v0, Vector3 &v1) {
    Vector3 v = v0 + 0.0042723144395936940576063989283284 * v1;
    integrand.weighPoint( 0.022935322010529224963732008058970,v);
    v = v0 + 0.025446043828620756865888097308925 * v1;
    integrand.weighPoint(  0.063092092629978553290700663189204,v);
    v = v0 + 0.067567788320115451661251881887438 * v1;
    integrand.weighPoint(  0.104790010322250183839876322541518,v);
    v = v0 + 0.12923440720030276995800022632466 * v1;
    integrand.weighPoint(  0.140653259715525918745189590510238,v);
    v = v0 + 0.20695638226615442611944217787823 * v1;
    integrand.weighPoint(  0.169004726639267902826583426598550,v);
    v = v0 + 0.29707742431130140792205907018797 * v1;
    integrand.weighPoint(  0.190350578064785409913256402421014,v);
    v = v0 + 0.3961075224960507457083735971537 * v1;
    integrand.weighPoint(  0.204432940075298892414161999234649,v);
    v = v0 + 0.5 * v1;
    integrand.weighPoint(  0.209482141084727828012999174891714,v);
    v = v0 + 0.6038924775039492542916264028463 * v1;
    integrand.weighPoint(  0.204432940075298892414161999234649,v);
    v = v0 + 0.7029225756886985365667896985542 * v1;
    integrand.weighPoint(  0.190350578064785409913256402421014,v);
    v = v0 +  0.79304361773384557388055782212177 * v1;
    integrand.weighPoint(  0.169004726639267902826583426598550,v);
    v = v0 +   0.87076559279969723004199977367534 * v1;
    integrand.weighPoint(  0.140653259715525918745189590510238,v);
    v = v0 +   0.93243221167988454833874811811256 * v1;
    integrand.weighPoint(  0.104790010322250183839876322541518,v);
    v = v0 +   0.97455395617137918762296067143325 * v1;
    integrand.weighPoint(  0.063092092629978553290700663189204,v);
    v = v0 +   0.99572768556040625043124236981384 * v1;
    integrand.weighPoint(  0.022935322010529224963732008058970,v);
  }
};




#endif
