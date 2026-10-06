#pragma once

#include "Common/Float.h"
#include "AR/variant.h"

namespace AR {

    //! force class of one particle for LogH
    /*! Store one particle's acceleration and perturbation
     */
    struct ForceLogH {
    public:
        Float acc_in[3];   ///< total acceleration from all inner particle members of the group
        Float acc_pert[3]; ///< perturbation from outside of the group
        Float pot_in;      ///< inner potential
        Float pot_pert;    ///< perturbation potential

        //! initialization 
        /*! set all data to zero
         */
        ForceLogH(): acc_in{Float(0.0),Float(0.0),Float(0.0)}, 
                     acc_pert{Float(0.0),Float(0.0),Float(0.0)},
                     pot_in(0.0), pot_pert(0.0) {}

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            fwrite(this, sizeof(*this),1,_fout);
        }

        //! read class data with BINARY format 
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this), 1, _fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }

        //! clear
        void clear() {
            acc_in[0] = acc_in[1] = acc_in[2] = 0.0;
            acc_pert[0] = acc_pert[1] = acc_pert[2] = 0.0;
            pot_in = pot_pert = 0.0;
        }
    };

    //! force class of one particle for TTL
    /*! Store one particle's acceleration, perturbation and time transformation function gradient
     */
    struct ForceTTL {
    public:
        Float acc_in[3];   ///< total acceleration from all inner particle members of the group
        Float acc_pert[3]; ///< perturbation from outside of the group
        Float pot_in;      ///< inner potential
        Float pot_pert;    ///< perturbation potential
        Float gtgrad[3];   ///< time transformation function gradient

        //! initialization 
        /*! set all data to zero
         */
        ForceTTL(): acc_in{Float(0.0),Float(0.0),Float(0.0)}, 
                    acc_pert{Float(0.0),Float(0.0),Float(0.0)},
                    pot_in(0.0), pot_pert(0.0),
                    gtgrad{Float(0.0),Float(0.0),Float(0.0)} {}

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            fwrite(this, sizeof(*this),1,_fout);
        }

        //! read class data with BINARY format 
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this), 1, _fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }

        //! clear
        void clear() {
            acc_in[0] = acc_in[1] = acc_in[2] = 0.0;
            acc_pert[0] = acc_pert[1] = acc_pert[2] = 0.0;
            pot_in = pot_pert = 0.0;
            gtgrad[0] = gtgrad[1] = gtgrad[2] = 0.0;
        }
    };

    //! the force type is determined by the time transformation: TTL adds the gradient of the time transformation function
    template <class Ttransform> struct ForceTypeOf;
    template <> struct ForceTypeOf<LogH> { typedef ForceLogH type; };
    template <> struct ForceTypeOf<TTL> { typedef ForceTTL type; };

    template <class Ttransform>
    using ForceType = typename ForceTypeOf<Ttransform>::type;

    //! force class of one particle for the default time transformation
    typedef ForceType<DefaultTimeTransformation> Force;

}
