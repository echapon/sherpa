#include "AMEGIC++/Main/Point.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace MODEL;
using namespace ATOOLS;

Point::Point(const Point& copy) { 
  Color   = NULL;
  Lorentz = NULL;
  middle  = 0;

  *this = copy;
} 

Point::Point(int extra) { 
  zwf     = 0;
  propid  = 0;
  v       = NULL;
  Color   = NULL;
  Lorentz = NULL;
  middle  = 0;
}

Point& Point::operator=(const Point& p) {
  if (this!=&p) {
    number = p.number;
    b      = p.b;
    t      = p.t;
    zwf    = p.zwf;
    propid = p.propid;
    m      = p.m;
    fl     = p.fl;

    if (p.Lorentz) {
      if (Color==NULL) Color = new Color_Function();
      *Color = *p.Color; 
      if (Lorentz) Lorentz->Delete();
      Lorentz = p.Lorentz->GetCopy(); 
    }

    left   = p.left;
    right  = p.right;
    middle = p.middle;
    prev  = p.prev;
    v = p.v;
    //cpl's
    cpl.clear();

    for(size_t i=0;i<p.Ncpl();i++) cpl.push_back(p.cpl[i]);
  }
  return *this;
}

void Point::ResetExternalNumbers(int os)
{
  if (number<100 && b==1) number+=os;
  if (left) {
    left->ResetExternalNumbers(os);
    right->ResetExternalNumbers(os);
    if (middle) middle->ResetExternalNumbers(os);
  }
}

void Point::ResetFlag()
{
  t = 0;
  if (left) {
    left->ResetFlag();
    right->ResetFlag();
    if (middle) middle->ResetFlag();
  }
}

void Point::ResetProps()
{
  int st = 0;
  ResetProps(st);
}

void Point::ResetProps(int &st)
{
  if (b==2) b=1;
  if (left) {
    if (number!=0){
      st++;
      number = st;
      if (fl.IsFermion()) number+=100;
      if (fl.IsBoson())   number+=200;
    }
    left->ResetProps(st);
    right->ResetProps(st);
    if (middle) middle->ResetProps(st);
  }
}

Point* Point::CopyList(Point* p)
{
  *this = p[0];
  Point* nx = this;
  if (p[0].left) {
    left = nx + 1;
    right = left->CopyList(p[0].left) + 1;
    nx = right->CopyList(p[0].right);
    if (p[0].middle) {
      middle = nx + 1;
      nx = middle->CopyList(p[0].middle);
    }
  }
  return nx;
}

int Point::CountKK()
{
  int KKnum=0;
  if (left) {
    KKnum+=left->CountKK();
    KKnum+=right->CountKK();
    if (middle) KKnum+=middle->CountKK();
  }
  if (fl.IsKK()) KKnum++;
  return KKnum;
}

bool Point::CountT(int & tchan,const long unsigned int & kfcode) {
  long unsigned int comp;
  if (left) {
    if (left->CountT(tchan,kfcode)) { 
      comp=left->fl.Kfcode();
      if ((comp==kfcode || kfcode==0) && !fl.Strong()) tchan++; 
      return true;
    }
    if (right->CountT(tchan,kfcode)) {
      comp=right->fl.Kfcode();
      if ((comp==kfcode || kfcode==0) && !fl.Strong()) tchan++; 
      return true;
    }
    if (middle && middle->CountT(tchan,kfcode)) {
      comp=middle->fl.Kfcode();
      if ((comp==kfcode || kfcode==0) && !fl.Strong()) tchan++; 
      return true;
    }
  }
  else if (b==-1) return true;
  return false;
}

void Point::GeneratePropID()
{
  propid=0;
  if (!left) {
    propid=(1<<number);
    return;
  }
  left->GeneratePropID();
  propid+=left->propid;
  right->GeneratePropID();
  propid+=right->propid;
  if (middle) {
    middle->GeneratePropID();
    propid+=middle->propid;    
  }
}

std::string Point::GetPropID() const
{
  return fl.IDName()+ATOOLS::ToString(propid);
}

void Point::FindOrder(std::vector<int> &order)
{
  // HS added this to prevent run-time segfault when compiled with gcc6
  if (v == NULL) return;

  if (v) {
    if (order.size()<v->order.size())        order.resize(v->order.size(),0);
    for (size_t i(0);i<v->order.size();++i)  order[i]+=v->order[i];
  }
  if (left) left->FindOrder(order);
  if (right) right->FindOrder(order);
  if (middle) middle->FindOrder(order);
}

std::ostream &AMEGIC::operator<<(std::ostream &str,const Point &p)
{
  str<<p.fl<<"("<<p.b<<","<<p.number;
  if (p.v) str<<",order="<<p.v->order;
  if (p.Color) str<<",col="<<*p.Color;
  if (p.Lorentz) str<<",lorentz="<<*p.Lorentz;
  str<<")";
  if (p.left) {
    str<<"[->"<<*p.left<<","<<*p.right;
    if (p.middle) str<<","<<*p.middle;
    str<<"]";
  }
  return str;
}

