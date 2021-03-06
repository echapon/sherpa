#ifndef ATOOLS_Phys_Blob_List_H
#define ATOOLS_Phys_Blob_List_H

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Event_Weights.H"
#include "ATOOLS/Phys/Particle_List.H"
#include <deque>
#include <map>
#include <set>
#include <iostream>

namespace ATOOLS {

  class Blob_List: public std::deque<Blob*> {
  private:
    static std::map<btp::code, unsigned long int> s_momfails;

    class Destructor {
    private:
      Blob_List *p_list;

      void RegisterDefaults() const;

    public:

      // destructor
      Destructor(): p_list(NULL) {}
      Destructor(Blob_List *const list): p_list(list) {}

      // destructor
      inline ~Destructor() { if (p_list!=NULL) p_list->Clear(); }

    };// end of class Destructor

    double m_extweight;
    Destructor m_destructor;

    friend std::ostream &operator<<(std::ostream &str,const Blob_List &list);

    void FindConnected(Blob *blob,Blob_List &connected,
		       std::set<const Blob*> &selected);

    void DeleteConnected(Blob *blob,std::set<Blob*> &deleted);

    bool TotalFourMomentum(Blob *blob,std::set<Blob*> &summed,
			   Vec4D &inisum,Vec4D &finsum,const int mode) const;

  public:

    // constructor
    Blob_List();
    Blob_List(const bool destruct);

    // member functions
    Vec4D TotalFourMomentum() const;
    Vec4D IncomingFourMomentum() const;
    Vec4D OutgoingFourMomentum() const;

    bool FourMomentumConservation() const;
    bool ColorConservation() const;

    Blob *FindFirst(const btp::code code) const;
    Blob *FindLast(const btp::code code) const;

    Blob_List Find(const btp::code code) const;

    Blob_List FindConnected(const Blob *blob);
    Blob_List FindConnected(const Particle *blob);

    Blob *AddBlob(const btp::code &type);

    bool Delete(Blob *blob);

    size_t DeleteConnected(Blob *blob);
    size_t DeleteConnected(Particle *particle);

    Particle_List ExtractParticles(const int status,
				   const int mode=0) const;
    Particle_List ExtractLooseParticles(const int mode=0) const;

    bool MergeSubsequentType(btp::code,btp::code,long int &,long int &);
    void MergeSubsequentTypeRecursively(btp::code,btp::code,
					long int &,long int &);
    void Clear(Blob *blob=NULL);

    Blob_List Copy() const;

    ATOOLS::Event_Weights Weights() const;
    double Weight() const;
    void   SetExternalWeight(const double & wt) { m_extweight = wt; }
    double ExternalWeight() const { return m_extweight; }
    static void PrintMomFailStatistics(std::ostream &str);
  };// end of class Blob_List
  
  std::ostream &operator<<(std::ostream &str,const Blob_List &list);

}// end of namespace ATOOLS

#endif
