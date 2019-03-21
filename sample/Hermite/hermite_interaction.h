#pragma once

#include "Common/Float.h"

//! hermite interaction class 
class HermiteInteraction{
public:
    Float eps_sq; // softening parameter
    Float G;      // gravitational constant

    // constructor
    HermiteInteraction(): eps_sq(Float(-1.0)), G(Float(-1.0)) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(eps_sq>=0.0);
        ASSERT(G>0.0);
        return true;
    }        

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"eps_sq: "<<eps_sq<<std::endl
             <<"G     : "<<G<<std::endl;
    }    

    //! calculate separation square between i and j particles
    /*!
      @param[in]: _pi: particle i
      @param[in]: _pj: particle j
      \return the distance square of the pair
     */
    template<class Tpi, class Tpj>
    inline Float calcR2Pair(const Tpi& _pi,
                            const Tpj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        return dr2;
    }

    //! calculate acceleration and jerk of one pair
    /*!
      @param[out]: _fi: acceleration for i particle
      @param[in]: _pi: particle i
      @param[in]: _pj: particle j
      \return the distance square of the pair
     */
    template<class Tpi, class Tpj>
    inline Float calcAccJerkPairSingleSingle(H4::ForceH4& _fi,
                                             const Tpi& _pi,
                                             const Tpj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float dv[3] = {_pj.vel[0] - _pi.vel[0],
                             _pj.vel[1] - _pi.vel[1],
                             _pj.vel[2] - _pi.vel[2]};
        const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
        const Float r = sqrt(dr2_eps);
        const Float rinv = 1.0/r;
            
        const Float rinv2 = rinv*rinv;
        const Float rinv3 = rinv2*rinv;

        const Float mor3 = G*_pj.mass*rinv3;
        const Float acc0[3] = {mor3*dr[0], mor3*dr[1], mor3*dr[2]};
        const Float acc1[3] = {mor3*dv[0] - 3.0*drdv*rinv2*acc0[0],
                               mor3*dv[1] - 3.0*drdv*rinv2*acc0[1],
                               mor3*dv[2] - 3.0*drdv*rinv2*acc0[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        return dr2;
    }

    //! calculate acceleration and jerk of one pair single and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      \return the minimum distance square of i and member j
     */
    template<class Tpi, class Tgroup>
    inline Float calcAccJerkPairSingleGroupMember(H4::ForceH4& _fi,
                                                  const Tpi& _pi,
                                                  const Tgroup& _gj) {
        const int n_member = _gj.particles.getSize();
        auto* member_adr = _gj.particles.getOriginAddressArray();
        Float r2_min = NUMERIC_FLOAT_MAX;
        for (int k=0; k<n_member; k++) {
            const auto& pj = *member_adr[k];
            ASSERT(_pi.id!=pj.id);
            const Float dr[3] = {pj.pos[0]-_pi.pos[0], 
                                 pj.pos[1]-_pi.pos[1],
                                 pj.pos[2]-_pi.pos[2]};
            Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            Float dr2_eps = dr2 + eps_sq;
            const Float dv[3] = {pj.vel[0] - _pi.vel[0],
                                 pj.vel[1] - _pi.vel[1],
                                 pj.vel[2] - _pi.vel[2]};
            const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
            const Float r = sqrt(dr2_eps);
            const Float rinv = 1.0/r;
            
            const Float rinv2 = rinv*rinv;
            const Float rinv3 = rinv2*rinv;

            const Float mor3 = G*pj.mass*rinv3;
            const Float acc0[3] = {mor3*dr[0], mor3*dr[1], mor3*dr[2]};
            const Float acc1[3] = {mor3*dv[0] - 3.0*drdv*rinv2*acc0[0],
                                   mor3*dv[1] - 3.0*drdv*rinv2*acc0[1],
                                   mor3*dv[2] - 3.0*drdv*rinv2*acc0[2]};
            _fi.acc0[0] += acc0[0];
            _fi.acc0[1] += acc0[1];
            _fi.acc0[2] += acc0[2];

            _fi.acc1[0] += acc1[0];
            _fi.acc1[1] += acc1[1];
            _fi.acc1[2] += acc1[2];

            r2_min = std::min(r2_min, dr2);
        }
        return r2_min;
    }

    //! calculate acceleration and jerk of one pair single and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      @param[in]: _pj: predicted group j cm
      \return the distance square of the pair
     */
    template<class Tpi, class Tgroup, class Tpcmj>
    inline Float calcAccJerkPairSingleGroupCM(H4::ForceH4& _fi,
                                              const Tpi& _pi,
                                              const Tgroup& _gj,
                                              const Tpcmj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float dv[3] = {_pj.vel[0] - _pi.vel[0],
                             _pj.vel[1] - _pi.vel[1],
                             _pj.vel[2] - _pi.vel[2]};
        const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
        const Float r = sqrt(dr2_eps);
        const Float rinv = 1.0/r;
            
        const Float rinv2 = rinv*rinv;
        const Float rinv3 = rinv2*rinv;

        const Float mor3 = G*_pj.mass*rinv3;
        const Float acc0[3] = {mor3*dr[0], mor3*dr[1], mor3*dr[2]};
        const Float acc1[3] = {mor3*dv[0] - 3.0*drdv*rinv2*acc0[0],
                               mor3*dv[1] - 3.0*drdv*rinv2*acc0[1],
                               mor3*dv[2] - 3.0*drdv*rinv2*acc0[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        return dr2;
    }

    //! calculate acceleration and jerk of one pair group cm and single
    /*!
      @param[out]: _fi: acceleration for i particle
      @param[in]: _gi: particle group i
      @param[in]: _pi: particle i
      @param[in]: _pj: particle j
      \return the distance square of the pair
     */
    template<class Tgroup, class Tpcmi, class Tpj>
    inline Float calcAccJerkPairGroupCMSingle(H4::ForceH4& _fi,
                                              const Tgroup& _gi,
                                              const Tpcmi& _pi,
                                              const Tpj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float dv[3] = {_pj.vel[0] - _pi.vel[0],
                             _pj.vel[1] - _pi.vel[1],
                             _pj.vel[2] - _pi.vel[2]};
        const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
        const Float r = sqrt(dr2_eps);
        const Float rinv = 1.0/r;
            
        const Float rinv2 = rinv*rinv;
        const Float rinv3 = rinv2*rinv;

        const Float mor3 = G*_pj.mass*rinv3;
        const Float acc0[3] = {mor3*dr[0], mor3*dr[1], mor3*dr[2]};
        const Float acc1[3] = {mor3*dv[0] - 3.0*drdv*rinv2*acc0[0],
                               mor3*dv[1] - 3.0*drdv*rinv2*acc0[1],
                               mor3*dv[2] - 3.0*drdv*rinv2*acc0[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        return dr2;
    }

    //! calculate acceleration and jerk of one pair group cm and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _gi: particle group i
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      \return the minimum distance square of i and member j
     */
    template<class Tpi, class Tgroup>
    inline Float calcAccJerkPairGroupCMGroupMember(H4::ForceH4& _fi,
                                                   const Tgroup& _gi,
                                                   const Tpi& _pi,
                                                   const Tgroup& _gj) {
        const int n_member = _gj.particles.getSize();
        auto* member_adr = _gj.particles.getOriginAddressArray();
        Float r2_min = NUMERIC_FLOAT_MAX;
        for (int k=0; k<n_member; k++) {
            const auto& pj = *member_adr[k];
            ASSERT(_pi.id!=pj.id);
            const Float dr[3] = {pj.pos[0]-_pi.pos[0], 
                                 pj.pos[1]-_pi.pos[1],
                                 pj.pos[2]-_pi.pos[2]};
            Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            Float dr2_eps = dr2 + eps_sq;
            const Float dv[3] = {pj.vel[0] - _pi.vel[0],
                                 pj.vel[1] - _pi.vel[1],
                                 pj.vel[2] - _pi.vel[2]};
            const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
            const Float r = sqrt(dr2_eps);
            const Float rinv = 1.0/r;
            
            const Float rinv2 = rinv*rinv;
            const Float rinv3 = rinv2*rinv;

            const Float mor3 = G*pj.mass*rinv3;
            const Float acc0[3] = {mor3*dr[0], mor3*dr[1], mor3*dr[2]};
            const Float acc1[3] = {mor3*dv[0] - 3.0*drdv*rinv2*acc0[0],
                                   mor3*dv[1] - 3.0*drdv*rinv2*acc0[1],
                                   mor3*dv[2] - 3.0*drdv*rinv2*acc0[2]};
            _fi.acc0[0] += acc0[0];
            _fi.acc0[1] += acc0[1];
            _fi.acc0[2] += acc0[2];

            _fi.acc1[0] += acc1[0];
            _fi.acc1[1] += acc1[1];
            _fi.acc1[2] += acc1[2];

            r2_min = std::min(r2_min, dr2);
        }
        return r2_min;
    }

    //! calculate acceleration and jerk of one pair single and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _gi: particle group i
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      @param[in]: _pj: predicted group j cm
      \return the distance square of the pair
     */
    template<class Tpcmi, class Tgroup, class Tpcmj>
    inline Float calcAccJerkPairGroupCMGroupCM(H4::ForceH4& _fi,
                                               const Tgroup& _gi,
                                               const Tpcmi& _pi,
                                               const Tgroup& _gj,
                                               const Tpcmj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float dv[3] = {_pj.vel[0] - _pi.vel[0],
                             _pj.vel[1] - _pi.vel[1],
                             _pj.vel[2] - _pi.vel[2]};
        const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
        const Float r = sqrt(dr2_eps);
        const Float rinv = 1.0/r;
            
        const Float rinv2 = rinv*rinv;
        const Float rinv3 = rinv2*rinv;

        const Float mor3 = G*_pj.mass*rinv3;
        const Float acc0[3] = {mor3*dr[0], mor3*dr[1], mor3*dr[2]};
        const Float acc1[3] = {mor3*dv[0] - 3.0*drdv*rinv2*acc0[0],
                               mor3*dv[1] - 3.0*drdv*rinv2*acc0[1],
                               mor3*dv[2] - 3.0*drdv*rinv2*acc0[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        return dr2;
    }

    //! calculate potential from j to i
    /*! \return the potential
     */
    template<class Tpi, class Tpj>
    inline Float calcPotPair(const Tpi& _pi,
                             const Tpj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float r = sqrt(dr2_eps);
        const Float rinv = 1.0/r;
        
        return -G*_pj.mass*rinv;
    }
    
    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(*this),1,_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }    
};

