#ifndef LEVEL_H
#define LEVEL_H

#include <stdexcept>
#include <vector>
#include <unordered_set>
#include <array>
#include <unordered_map>
#include <algorithm>

#include "SphericalCoordinates.h"

#include <Eigen/Dense>

class BoxTree;

class Level {

    friend class BoxTree;

    private:

        const int level_;

        const double min_x_;
        const double min_y_;
        const double min_z_;

        const double boxsize_;

        const double wavenumber_;

        const long long nboxesperside_;

        double scaling_r2s_;
        long long nconeselevation_;
        std::vector<double> radial_intervals_cone_;

        std::unordered_set<long long> mortonidofrelboxes_;
        std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_;

        std::unordered_set<long long> mortonidofrelboxes_2_;
        std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_2_;

        std::vector<long long> relconesegment2cocenteredmortonboxid_;
        std::vector<long long> relconesegment2nonrelconesegment_;
        std::unordered_map<long long, std::unordered_map<long long, long long>> mortonboxnonrelcone2relcone_;  

        std::vector<long long> split_points_relconesegments_;  
        std::unordered_map<long long, std::unordered_map<long long, long long>> mortonboxnonrelcone2relcone_2_;

        std::unordered_set<long long> propagationrequiredrelconesegments_;
        std::unordered_set<long long> interpolationrequiredrelconesegments_;   

        std::vector<std::vector<long long>> propagation_not_in_rank_;
        std::vector<int> recv_counts_propagation_;
        std::vector<int> displs_propagation_;
        std::vector<long long> cone_segments_propagation_not_in_rank_all_;
        std::vector<std::vector<long long>> interpolation_not_in_rank_;
        std::vector<int> recv_counts_interpolation_;
        std::vector<int> displs_interpolation_;
        std::vector<long long> cone_segments_interpolation_not_in_rank_all_;

        std::unordered_map<long long, std::array<std::unordered_map<long long, std::vector<int>>, 8>> protobox_; 
          
        std::unordered_map<long long, std::array<double, 3>> means_;

        std::vector<double> relconesegment2minradius_;
        std::unordered_map<long long, double> mortonboxnonrelcone2minradius_;
        std::vector<double> relconesegment2maxradius_;
        std::unordered_map<long long, double> mortonboxnonrelcone2maxradius_;

        std::unordered_map<long long, Eigen::MatrixXd> interpolation_points_;
        std::vector<long long> indexes_relconesegment_;
        std::vector<long long> num_elems_indexes_relconesegment_;

    public:

        Level(const int level,
              const double min_x, const double min_y, const double min_z,
              const double boxsize, const std::complex<double> wavenumber) :
              level_(level), min_x_(min_x), min_y_(min_y), min_z_(min_z),
              boxsize_(boxsize), nboxesperside_(1 << level_), wavenumber_(wavenumber.real()) {

            InitializeConeSegments();           

        }

        ~Level() {

        }

        void InitializeConeSegments() {

            const double minsize = (3.0/2.0 * boxsize_) * (1.0-1e-10);
            const double maxsize = (std::sqrt(3.0) * boxsize_ * nboxesperside_ - std::sqrt(3.0)/2.0 * boxsize_) * (1.0+1e-10);
            constexpr long long nconesradial = 1;
            constexpr long long nconeselevation = 2; // elevation goes from 0 to pi
            long long nconeslevelscaling = 1;

            if (level_ > 1)
                nconeslevelscaling = std::max<long long>(1, std::ceil(wavenumber_ * boxsize_ * M_1_PI * 1.2)); //?

            nconeselevation_ = nconeslevelscaling * nconeselevation;
            const long long nconesradialonlevel = nconeslevelscaling * nconesradial;

            scaling_r2s_ = std::sqrt(3.0)/2.0 * boxsize_;
            const double eta = 1.0/std::sqrt(3.0)-1e-5;

            const double tmpsizeinr = (maxsize - minsize + 1e-10);
            const double tmpsizeins = std::min(eta/nconesradialonlevel, scaling_r2s_*(1.0/minsize - 1.0/(minsize+tmpsizeinr)));
            
            const long long nactualcones = std::ceil(scaling_r2s_*(1.0/minsize - 1.0/maxsize)/tmpsizeins);
            radial_intervals_cone_.resize(nactualcones+1);
            radial_intervals_cone_[0] = minsize;

            const double finalsizeins = (scaling_r2s_*(1.0/minsize - 1.0/maxsize)/nactualcones);
            for (long long i = nactualcones-1; i >= 0; i--) {
                const long long iter = nactualcones - i;       
                radial_intervals_cone_[iter] = scaling_r2s_*maxsize/(scaling_r2s_ + maxsize*i*finalsizeins);
            }

        }

        static long long Point2Morton (const double x, const double y, const double z,
                                       const double min_x, const double min_y, const double min_z,
                                       const double boxsize, const int level) {

            const long long posx = static_cast<long long>((x - min_x)/boxsize);
            const long long posy = static_cast<long long>((y - min_y)/boxsize);
            const long long posz = static_cast<long long>((z - min_z)/boxsize);

            if (posx < 0 || posy < 0 || posz < 0 || posx >= 1 << level || posy >= 1 << level || posz >= 1 << level)
                throw std::logic_error("The 3D box index cannot be < 0 or larger than the number of boxes");

            return Box2Morton(posx, posy, posz, level);

        }

        long long  Point2Morton(const double x, const double y, const double z) const {

            return Point2Morton(x, y, z, min_x_, min_y_, min_z_, boxsize_, level_);

        }

        static long long Box2Morton(const long long x, const long long y, const long long z, const int level) {

            uint64_t answer = 0;
            answer |= splitBy3(x, level) | splitBy3(y, level) << 1 | splitBy3(z, level) << 2;

            return answer;

        }

        long long Box2Morton(const long long x, const long long y, const long long z) const {

            return Box2Morton(x, y, z, level_);

        }

        static uint64_t splitBy3(const unsigned int a, const int level) {

            uint64_t x = a & 0x1fffff;
            x &= (static_cast<uint64_t>(1) << (level+1)) -1;
            x = (x | x << 32) & 0x1f00000000ffff; 
            x = (x | x << 16) & 0x1f0000ff0000ff;
            x = (x | x << 8) & 0x100f00f00f00f00f;
            x = (x | x << 4) & 0x10c30c30c30c30c3;
            x = (x | x << 2) & 0x1249249249249249;

            return x;

        }

        inline uint64_t splitBy3(const unsigned int a) const {

            return splitBy3(a, level_);

        }

        static unsigned int mergeFrom3(uint64_t x) {

            x &= 0x1249249249249249;
            x = (x | x >> 2) & 0x10c30c30c30c30c3;
            x = (x | x >> 4) & 0x100f00f00f00f00f;
            x = (x | x >> 8) & 0x1f0000ff0000ff;
            x = (x | x >> 16) & 0x1f00000000ffff;
            x = (x | x >> 32) & 0x1fffff;

            return x;

        }

        static void Morton2Box(long long morton_box, long long& x, long long& y, long long& z, const int level) {

            if (morton_box >= std::pow(2.0, 3*level)) {
                std::cout << morton_box << " " << " " << level << " " << (1 >> (3*level)) << " " << std::pow(2.0, 3*level) << "\n";
                throw std::invalid_argument("Current level cannot have this morton order.");
            }

            x = mergeFrom3(morton_box);
            y = mergeFrom3(morton_box >> 1);
            z = mergeFrom3(morton_box >> 2);

        }

        void Morton2Box(long long morton_box, long long& x, long long& y, long long& z) const {

            Morton2Box(morton_box, x, y, z, level_);

        }

        static long long even (const long long x) {

            return static_cast<long long>(x/2)*2;

        }

        bool IsRelevant(long long morton_box) const {

            return mortonidofrelboxes_2_.count(morton_box) > 0;

        }

        std::vector<long long> GetCousins(const long long morton_box) const {

            std::vector<long long> cousins;
            cousins.reserve(189);

            long long x, y, z;
            Morton2Box(morton_box, x, y, z);
            long long xlo = std::max<long long>(even(x)-2, 0);
            long long xup = std::min<long long>(nboxesperside_-1, even(x)+3);
            long long ylo = std::max<long long>(even(y)-2, 0);
            long long yup = std::min<long long>(nboxesperside_-1, even(y)+3);
            long long zlo = std::max<long long>(even(z)-2, 0);
            long long zup = std::min<long long>(nboxesperside_-1, even(z)+3);
            
            for (int k = zlo; k <= zup; k++) {
                for (int j = ylo; j <= yup; j++) {
                    for (int i = xlo; i <= xup; i++) {
                        if (std::abs(i - x) > 1 || 
                            std::abs(j - y) > 1 ||
                            std::abs(k - z) > 1) {

                                const long long cousin_morton = Box2Morton(i, j, k);
                                if (IsRelevant(cousin_morton))
                                    cousins.push_back(cousin_morton);

                        }
                    }
                }
            }

            std::sort(cousins.begin(), cousins.end());
            return cousins;

        }

        void GetBoxCenter(long long morton_box, double& x, double& y, double& z) const {

            long long posx, posy, posz;
            Morton2Box(morton_box, posx, posy, posz);
            x = min_x_ + posx * boxsize_ + 0.5 * boxsize_;
            y = min_y_ + posy * boxsize_ + 0.5 * boxsize_;
            z = min_z_ + posz * boxsize_ + 0.5 * boxsize_;

        }

        long long GetNRadialIntervals() const  {return radial_intervals_cone_.size() - 1;}

        long long Threedimensionalconesegment2nonrelconesegment(long long r, long long theta, long long phi) const {

            return phi * nconeselevation_ * GetNRadialIntervals() + theta * GetNRadialIntervals() + r;

        }

        long long Cart2Nonrelcone(const long long mortonbox, double x, double y, double z) const {

            double centerx, centery, centerz;
            GetBoxCenter(mortonbox, centerx, centery, centerz);
            x -= centerx;
            y -= centery;
            z -= centerz;
            Functions::CartToSph(x, y, z);
            const long long posr = static_cast<long long>(std::upper_bound(radial_intervals_cone_.begin(), radial_intervals_cone_.end(), x)
            - radial_intervals_cone_.begin()) - 1;
            const long long postheta = static_cast<long long>(y * nconeselevation_ / M_PI);
            const long long posphi = static_cast<long long>(z * nconeselevation_ / M_PI);
            return Threedimensionalconesegment2nonrelconesegment(posr, postheta, posphi);

        }

        long long MortonboxNonrelconesegment2Relconesegment(long long morton_box, long long nonrelconesegment) const {

            return mortonboxnonrelcone2relcone_2_.at(morton_box).at(nonrelconesegment);
            
        }

        long long Nonrelconesegment2radialpos(const long long nonrelconesegment) const {

            return nonrelconesegment % GetNRadialIntervals();

        }

        void Nonrelconesegment2threedimensionalconesegment(const long long nonrelconesegment, long long & r, long long & theta, long long & phi) const {

            r = Nonrelconesegment2radialpos(nonrelconesegment);

            theta = nonrelconesegment / GetNRadialIntervals() % nconeselevation_;

            phi = nonrelconesegment / (nconeselevation_ * GetNRadialIntervals());

        }

        double ChebPosr2Radius(long long nonrelposinrdirection, const double r) const {

            double min = radial_intervals_cone_[nonrelposinrdirection];

            double reciprocal_min = scaling_r2s_ / radial_intervals_cone_[nonrelposinrdirection+1];
            double dr = scaling_r2s_/min-reciprocal_min;

            return scaling_r2s_/((r + 1.0)/2.0*dr + reciprocal_min);

        }

        void Cheb2Cart(const long long mortonbox, long long nonrelconesegment, double & x, double & y, double & z) const {

            long long posr, postheta, posphi;
            
            Nonrelconesegment2threedimensionalconesegment(nonrelconesegment, posr, postheta, posphi);

            const double dangle = M_PI/nconeselevation_;

            x = ChebPosr2Radius(posr, x);
            y = (y + 1.0)/2.0*dangle + dangle * postheta;
            z = (z + 1.0)/2.0*dangle + dangle * posphi;

            double centerx, centery, centerz;
            GetBoxCenter(mortonbox, centerx, centery, centerz);
            Functions::SphToCart(x, y, z);
            x += centerx;
            y += centery;
            z += centerz;

        }

        std::vector<long long> GetNeighbours(long long morton_box) {

            std::vector<long long> neighbours;
            neighbours.reserve(27);
            long long x, y, z;
            Morton2Box(morton_box, x, y, z);
            long long xlo = std::max<long long>(0, x-1);
            long long xup = std::min<long long>(nboxesperside_-1, x+1);
            long long ylo = std::max<long long>(0, y-1);
            long long yup = std::min<long long>(nboxesperside_-1, y+1);
            long long zlo = std::max<long long>(0, z-1);
            long long zup = std::min<long long>(nboxesperside_-1, z+1);

            for (int k = zlo; k <= zup; k++) {
                for (int j = ylo; j <= yup; j++) {
                    for (int i = xlo; i <= xup; i++) {

                        const long long neighbour_morton = Box2Morton(i, j, k);

                        if (IsRelevant(neighbour_morton))
                            neighbours.push_back(neighbour_morton);

                    }
                }
            }

            std::sort(neighbours.begin(), neighbours.end());
            return neighbours;

        }

        double Cheb2Radius(long long nonrelconesegment, const double r) const {

            const long long nonrelposindirection = Nonrelconesegment2radialpos(nonrelconesegment);
            return ChebPosr2Radius(nonrelposindirection, r);

        }

        double SphRadius2ChebRadius(long long nonrelconesegment_rpos, const double r) const {

            double min = radial_intervals_cone_[nonrelconesegment_rpos];

            double reciprocal_min = scaling_r2s_/radial_intervals_cone_[nonrelconesegment_rpos+1];
            double dr = scaling_r2s_/min-reciprocal_min;

            return (scaling_r2s_/r - reciprocal_min)/dr*2.0 - 1.0;
            
        }

        void Cart2Cheb(const long long mortonbox, const long long nonrelconesegment, double& x, double& y, double& z) const {

            double centerx, centery, centerz;
            GetBoxCenter(mortonbox, centerx, centery, centerz);
            x -= centerx;
            y -= centery;
            z -= centerz;
            Functions::CartToSph(x, y, z);
            long long posr, postheta, posphi;
            Nonrelconesegment2threedimensionalconesegment(nonrelconesegment, posr, postheta, posphi);
            const double dangle = M_PI/nconeselevation_;
            x = SphRadius2ChebRadius(posr, x);
            y = (y - dangle*postheta)/dangle*2.0 - 1.0;
            z = (z - dangle*posphi)/dangle*2.0 - 1.0;

        }

        void Cart2Cheb(const long long mortonbox, double& x, double& y, double& z, long long& nonrelconesegment_new) const {

            double centerx, centery, centerz;
            GetBoxCenter(mortonbox, centerx, centery, centerz);
            x -= centerx;
            y -= centery;
            z -= centerz;
            Functions::CartToSph(x, y, z);
            const long long posr = static_cast<long long>(std::upper_bound(radial_intervals_cone_.begin(), radial_intervals_cone_.end(), x) - radial_intervals_cone_.begin()) - 1;
            const long long postheta = static_cast<long long>(y * nconeselevation_ / M_PI);
            const long long posphi = static_cast<long long>(z * nconeselevation_ / M_PI);
            nonrelconesegment_new = Threedimensionalconesegment2nonrelconesegment(posr, postheta, posphi);
            const double dangle = M_PI/nconeselevation_;
            x = SphRadius2ChebRadius(posr, x);
            y = (y - dangle*postheta)/dangle*2.0 - 1.0;
            z = (z - dangle*posphi)/dangle*2.0 - 1.0;

        }

        std::vector<long long> GetAllNeighboursAndCousins(const long long morton_box) const
        {
            
            std::vector<long long> neighboursandcousins;

            neighboursandcousins.reserve(216);

            long long x, y, z;
            Morton2Box(morton_box, x, y, z);

            long long xlo = std::max<long long>(even(x)-2, 0);
            long long xup = std::min<long long>(nboxesperside_-1, even(x)+3);
            long long ylo = std::max<long long>(even(y)-2, 0);
            long long yup = std::min<long long>(nboxesperside_-1, even(y)+3);
            long long zlo = std::max<long long>(even(z)-2, 0);
            long long zup = std::min<long long>(nboxesperside_-1, even(z)+3);

            for (int k = zlo; k <= zup; k++) {
                for (int j = ylo; j <= yup; j++) {
                    for (int i = xlo; i <= xup; i++) {

                            neighboursandcousins.push_back(Box2Morton(i, j, k));

                    }
                }
            }

            return neighboursandcousins;

        }

        std::vector<long long> GetNeighboursAll(long long morton_box) {

            std::vector<long long> neighbours;
            neighbours.reserve(27);
            long long x, y, z;
            Morton2Box(morton_box, x, y, z);
            long long xlo = std::max<long long>(0, x-1);
            long long xup = std::min<long long>(nboxesperside_-1, x+1);
            long long ylo = std::max<long long>(0, y-1);
            long long yup = std::min<long long>(nboxesperside_-1, y+1);
            long long zlo = std::max<long long>(0, z-1);
            long long zup = std::min<long long>(nboxesperside_-1, z+1);

            for (int k = zlo; k <= zup; k++) {
                for (int j = ylo; j <= yup; j++) {
                    for (int i = xlo; i <= xup; i++) {

                        const long long neighbour_morton = Box2Morton(i, j, k);

                        neighbours.push_back(neighbour_morton);

                    }
                }
            }

            std::sort(neighbours.begin(), neighbours.end());
            return neighbours;

        }

};

#endif