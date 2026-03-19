//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file TurbineDomain.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_TURBINEDOMAIN_T
#define TURBINECORE_TURBINEDOMAIN_T

#pragma once

#include <blockforest/BlockForest.h>
#include <core/math/AABB.h>
#include <core/math/Vector3.h>
#include <memory>

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/math/AABB.h"

namespace turbine_core {

    namespace domain {

        HOST_DEVICE_PREFIX FORCEINLINE real_t sqDistanceLineToPoint( const real_t& pt, const real_t& min, const real_t& max  )
        {
            if (pt < min)
                return (min - pt) * (min - pt);
            if (pt > max)
                return (pt - max) * (pt - max);
            if( pt == max)
                return real_t(1e-6);
            return real_t(0);
        }

        HOST_DEVICE_PREFIX FORCEINLINE real_t sqDistancePointToAABB( const Vector3<real_t> & pt, const math::AABB& aabb )
        {
            real_t sq = 0.0;

            sq += sqDistanceLineToPoint( pt[0], aabb.xMin(), aabb.xMax() );
            sq += sqDistanceLineToPoint( pt[1], aabb.yMin(), aabb.yMax() );
            sq += sqDistanceLineToPoint( pt[2], aabb.zMin(), aabb.zMax() );

            return sq;
        }

        inline real_t sqDistancePointToAABBPeriodic( Vector3<real_t>  pt,
                                                     const AABB & aabb,
                                                     const walberla::math::AABB & domain,
                                                     const Vector3<bool> & periodic )
        {
            auto size = domain.sizes() * real_t(0.5);
            auto d = pt - aabb.center();

            if (periodic[0] && (d[0] < -size[0])) pt[0] += domain.sizes()[0];
            if (periodic[0] && (d[0] > +size[0])) pt[0] -= domain.sizes()[0];

            if (periodic[1] && (d[1] < -size[1])) pt[1] += domain.sizes()[1];
            if (periodic[1] && (d[1] > +size[1])) pt[1] -= domain.sizes()[1];

            if (periodic[2] && (d[2] < -size[2])) pt[2] += domain.sizes()[2];
            if (periodic[2] && (d[2] > +size[2])) pt[2] -= domain.sizes()[2];

            real_t sq = 0.0;

            sq += sqDistanceLineToPoint( pt[0], aabb.xMin(), aabb.xMax() );
            sq += sqDistanceLineToPoint( pt[1], aabb.yMin(), aabb.yMax() );
            sq += sqDistanceLineToPoint( pt[2], aabb.zMin(), aabb.zMax() );

            return sq;
        }

        inline bool isInsideAABB( const Vector3<real_t> & pt,
                                  const real_t radius,
                                  const AABB & aabb)
        {
            if (!aabb.contains(pt[0], pt[1], pt[2])) return false;
            if ((pt[0] - aabb.xMin()) < radius) return false;
            if ((aabb.xMax() - pt[0]) < radius) return false;
            if ((pt[1] - aabb.yMin()) < radius) return false;
            if ((aabb.yMax() - pt[1]) < radius) return false;
            if ((pt[2] - aabb.zMin()) < radius) return false;
            if ((aabb.zMax() - pt[2]) < radius) return false;
            return true;
        }

        class TurbineDomain
        {
        public:
            explicit TurbineDomain(const std::shared_ptr<walberla::blockforest::BlockForest>& blockForest,
                                   const std::vector<AABB> & turbineAABBs);

            /**
             * @brief If the BlockForest is changed this function has to be called in order to
             * update all internal caches!
             *
             * Updates the local caches for local and neighbor AABBs.
             */
            void refresh();

            bool isContainedInProcessSubdomain(uint_t rank, const Vector3<real_t> & pt) const;
            bool isContainedInLocalSubdomain(const Vector3<real_t>& pt, const real_t& radius) const;
            /// Is the sphere defined by \p pt and \p radius completely inside the local subdomain?
            /// \attention Also take into account periodicity!
            /// \param pt center of the sphere
            /// \param radius radius of the sphere
            bool   isContainedInProcessSubdomain(const Vector3<real_t>& pt, const real_t& radius) const;
            int    findContainingProcessRank(const Vector3<real_t>& pt) const;
            void   periodicallyMapToDomain(Vector3<real_t>& pt) const;
            std::vector<uint_t> getNeighborProcesses() const;
            bool   intersectsWithProcessSubdomain(uint_t rank, const Vector3<real_t> & pt, const real_t& radius) const;
            bool   intersectsWithProcessSubdomain(uint_t rank, const AABB & aabb) const;
            void   correctPointPosition(Vector3<real_t>& pt) const;

            const math::AABB& getUnionOfLocalAABBs() const {return unionOfLocalAABBs_;}
            size_t getNumLocalAABBs() const {return localAABBs_.size();}
            size_t getNumNeighborSubdomains() const {return neighborSubdomains_.size();}
            size_t getNumNeighborProcesses() const {return neighborProcesses_.size();}
            const std::vector<math::AABB> & localAABBs() const { return localAABBs_; }
            const std::vector<math::AABB> neighborAABBs( const uint_t rank ) const {

                WALBERLA_ASSERT(std::is_sorted(neighborSubdomains_.begin(),
                                               neighborSubdomains_.end(),
                                               [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank;}))

                std::vector<math::AABB> neighborAABBs;

                size_t idx = 0;
                WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                while (neighborSubdomains_[idx].rank != int(rank)) {
                    ++idx;
                    WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                }
                while (neighborSubdomains_[idx].rank == int(rank)) {
                    neighborAABBs.emplace_back(neighborSubdomains_[idx].aabb);
                    ++idx;
                    if (idx >= neighborSubdomains_.size()) break;
                    WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                }

                return neighborAABBs;

            }

        private:
            bool isInsideGlobalDomain(const Vector3<real_t>& pt, const real_t& radius) const;

            std::shared_ptr<walberla::blockforest::BlockForest> blockForest_;

            struct Subdomain {
                Subdomain(const int r, const walberla::BlockID& id, const math::AABB& ab) : rank(r), blockID(id), aabb(ab) {}
                int rank;
                walberla::BlockID blockID;
                AABB aabb;
            };

            int ownRank_ = -1;
            Vector3<bool> periodic_;

            std::vector<AABB>       localAABBs_;
            AABB                    unionOfLocalAABBs_;
            std::vector<Subdomain>  neighborSubdomains_;
            std::vector<uint_t>     neighborProcesses_;

            const std::vector<AABB> &     turbineAABBs_;
        };


        /// \post neighborSubdomains_ is sorted by rank
        TurbineDomain::TurbineDomain(const std::shared_ptr<walberla::blockforest::BlockForest>& blockForest,
                                     const std::vector<AABB> & turbineAABBs)
                : blockForest_(blockForest), turbineAABBs_(turbineAABBs)
        {
            refresh();
        }

        /// \post neighborSubdomains_ is sorted by rank
        void TurbineDomain::refresh() {

            ownRank_ = walberla::mpi::MPIManager::instance()->rank();

            periodic_[0] = blockForest_->isPeriodic(0);
            periodic_[1] = blockForest_->isPeriodic(1);
            periodic_[2] = blockForest_->isPeriodic(2);

            localAABBs_.clear();
            neighborSubdomains_.clear();
            neighborProcesses_.clear();
            unionOfLocalAABBs_ = AABB(Vector3<real_t>(real_t(0)), Vector3<real_t>(real_t(0)));

            if (blockForest_->empty()) return;

            for (auto& iBlk : *blockForest_) {

                const auto * blk = dynamic_cast<walberla::blockforest::Block*>(&iBlk);
                WALBERLA_ASSERT_NOT_NULLPTR( blk )

                auto blkAABB = blk->getAABB();
                bool containsTurbine = false;

                for(auto & turbineAABB : turbineAABBs_) {
                    if(turbineAABB.intersects(blkAABB)) {
                        containsTurbine = true;
                    }
                }

                if(!containsTurbine)
                    continue;

                localAABBs_.emplace_back(blk->getAABB());
                for (uint_t nb = 0; nb < blk->getNeighborhoodSize(); ++nb) {

                    if (int(blk->getNeighborProcess(nb)) == ownRank_) continue;

                    //check if neighbor aabb is already present
                    const walberla::BlockID& nbBlkId = blk->getNeighborId(nb);
                    if (std::find_if(neighborSubdomains_.begin(),
                                     neighborSubdomains_.end(),
                                     [&nbBlkId](const auto& subdomain){return subdomain.blockID == nbBlkId;}) ==
                        neighborSubdomains_.end())
                    {
                        auto neighbourAABB = AABB(blk->getNeighborAABB(nb));
                        bool neighbourContainsLocalTurbine = false;

                        for(auto & turbineAABB : turbineAABBs_) {
                            if(turbineAABB.intersects(neighbourAABB)) {
                                if(turbineAABB.intersects(localAABBs_.back())) {
                                    neighbourContainsLocalTurbine = true;
                                }
                            }
                        }
                        if(neighbourContainsLocalTurbine)
                            neighborSubdomains_.emplace_back(int(blk->getNeighborProcess(nb)), nbBlkId, neighbourAABB);
                    }
                }
            }

            for( auto & localAABB : localAABBs_ ) {
                unionOfLocalAABBs_.merge( localAABB );
            }

            //sort by rank
            std::sort(neighborSubdomains_.begin(),
                      neighborSubdomains_.end(),
                      [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank; });

            //generate list of neighbor processes
            int prevRank = -1;
            for (auto& subdomain : neighborSubdomains_) {
                if ((prevRank != subdomain.rank) && (subdomain.rank != ownRank_)) {
                    neighborProcesses_.emplace_back(uint_t(subdomain.rank));
                    prevRank = subdomain.rank;
                }
            }
        }

        bool TurbineDomain::isContainedInProcessSubdomain(const uint_t rank, const Vector3<real_t> & pt) const {

            if (blockForest_->empty()) return false;

            if (uint_t(ownRank_) == rank) {
                // check if point is in local subdomain by checking all aabbs
                for (auto& aabb : localAABBs_) {
                    if (aabb.contains(pt)) return true;
                }
            } else {

                WALBERLA_ASSERT(std::is_sorted(neighborSubdomains_.begin(),
                                               neighborSubdomains_.end(),
                                               [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank;}))

                auto begin = std::find_if(neighborSubdomains_.begin(),
                                          neighborSubdomains_.end(),
                                          [&rank](const auto& subdomain){ return subdomain.rank == int(rank); });
                if (begin == neighborSubdomains_.end()) return false; //no information for rank available
                auto end = std::find_if(begin,
                                        neighborSubdomains_.end(),
                                        [&rank](const auto& subdomain){ return subdomain.rank != int(rank); });

                for (auto it = begin; it != end; ++it) {
                    if (it->aabb.contains(pt)) return true;
                }
            }
            return false;
        }

        bool TurbineDomain::isContainedInLocalSubdomain(const Vector3<real_t> & pt,
                                                        const real_t& radius) const {
            return std::any_of(localAABBs_.begin(),
                               localAABBs_.end(),
                               [&](auto& aabb)
                               {return isInsideAABB(pt, radius, aabb);});
        }

        bool TurbineDomain::isContainedInProcessSubdomain(const Vector3<real_t> & pt, const real_t& radius) const {
            if (blockForest_->empty()) return false;

            //completely contained in local aabb?
            for (auto& aabb : localAABBs_) {
                if (aabb.contains(pt)) {
                    if (isInsideAABB(pt, radius, aabb)) {
                        return true;
                    }

                    break;
                }
            }

            //intersects one of the neighboring subdomains?
            return std::none_of(neighborSubdomains_.begin(),
                                neighborSubdomains_.end(),
                                [&](const auto &subdomain) {
                                    return sqDistancePointToAABB(pt, subdomain.aabb) < radius * radius;
                                });
        }

        int TurbineDomain::findContainingProcessRank(const Vector3<real_t> & pt) const {
            if (blockForest_->empty()) return -1;

            if (isContainedInProcessSubdomain(uint_t(ownRank_), pt)) return ownRank_;
            for( uint_t rank : getNeighborProcesses() ) {
                if (isContainedInProcessSubdomain(rank, pt)) return int(rank);
            }
            return -1;
        }
        void TurbineDomain::periodicallyMapToDomain(Vector3<real_t> & pt) const {
            blockForest_->mapToPeriodicDomain(pt[0], pt[1], pt[2]);
        }

        std::vector<uint_t> TurbineDomain::getNeighborProcesses() const {
            return neighborProcesses_;
        }

        bool TurbineDomain::intersectsWithProcessSubdomain(const uint_t rank, const Vector3<real_t> & pt, const real_t& radius) const {
            if (blockForest_->empty()) return false;

            if (uint_t(ownRank_) == rank) {
                //=====================
                // LOCAL DOMAIN
                if (isInsideGlobalDomain(pt, radius)) {
                    for (auto& aabb : localAABBs_) {
                        if (sqDistancePointToAABB(pt, aabb) <= radius * radius) return true;
                    }
                } else {
                    for (auto& aabb : localAABBs_) {
                        if (sqDistancePointToAABBPeriodic(pt, aabb, blockForest_->getDomain(), periodic_) <= radius * radius) return true;
                    }
                }
            } else {
                //=====================
                // NEIGHBORING DOMAIN
                WALBERLA_ASSERT(std::is_sorted(neighborSubdomains_.begin(),
                                               neighborSubdomains_.end(),
                                               [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank;}))

                if (isInsideGlobalDomain(pt, radius)) {
                    size_t idx = 0;
                    WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                    while (neighborSubdomains_[idx].rank != int(rank)) {
                        ++idx;
                        WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                    }
                    while (neighborSubdomains_[idx].rank == int(rank)) {
                        if (sqDistancePointToAABB(pt, neighborSubdomains_[idx].aabb) <= radius * radius) return true;
                        ++idx;
                        if (idx >= neighborSubdomains_.size()) break;
                        WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                    }
                } else {
                    size_t idx = 0;
                    WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                    while (neighborSubdomains_[idx].rank != int(rank)) {
                        ++idx;
                        WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                    }
                    while (neighborSubdomains_[idx].rank == int(rank)) {
                        if (sqDistancePointToAABBPeriodic(pt, neighborSubdomains_[idx].aabb, blockForest_->getDomain(), periodic_) <= radius * radius) return true;
                        ++idx;
                        if (idx >= neighborSubdomains_.size()) break;
                        WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                    }
                }
            }

            return false;
        }

        bool TurbineDomain::intersectsWithProcessSubdomain(const uint_t rank, const AABB & aabb) const {
            if (blockForest_->empty()) return false;

            if (uint_t(ownRank_) == rank) {
                //=====================
                // LOCAL DOMAIN
                for(auto & localAABB : localAABBs_) {
                    if( localAABB.intersects(aabb) )
                        return true;
                }
            } else {
                //=====================
                // NEIGHBORING DOMAIN
                WALBERLA_ASSERT(std::is_sorted(neighborSubdomains_.begin(),
                                               neighborSubdomains_.end(),
                                               [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank;}))

                Vector3<real_t> translation = blockForest_->getDomain().center();
                for(uint_t d = 0; d < 3; ++d) {
                    translation[d] *= int(periodic_[d]);
                }

                size_t idx = 0;
                WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                while (neighborSubdomains_[idx].rank != int(rank)) {
                    ++idx;
                    WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                }

                while (neighborSubdomains_[idx].rank == int(rank)) {
                    if( neighborSubdomains_[idx].aabb.intersects(aabb) ) return true;
                    if( neighborSubdomains_[idx].aabb.getTranslated(translation).intersects(aabb) ) return true;
                    if( neighborSubdomains_[idx].aabb.getTranslated(-translation).intersects(aabb) ) return true;

                    ++idx;
                    if (idx >= neighborSubdomains_.size()) break;
                    WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size())
                }
            }

            return false;
        }

        void TurbineDomain::correctPointPosition(Vector3<real_t> & pt) const {
            const Vector3<real_t>  center = unionOfLocalAABBs_.center();
            const Vector3<real_t>  dis = pt - center;

            const auto& domain = blockForest_->getDomain();

            if (periodic_[0] && (-domain.xSize() * 0.5 > dis[0])) pt[0] += domain.xSize();
            if (periodic_[0] && (+domain.xSize() * 0.5 < dis[0])) pt[0] -= domain.xSize();

            if (periodic_[1] && (-domain.ySize() * 0.5 > dis[1])) pt[1] += domain.ySize();
            if (periodic_[1] && (+domain.ySize() * 0.5 < dis[1])) pt[1] -= domain.ySize();

            if (periodic_[2] && (-domain.zSize() * 0.5 > dis[2])) pt[2] += domain.zSize();
            if (periodic_[2] && (+domain.zSize() * 0.5 < dis[2])) pt[2] -= domain.zSize();
        }


        bool TurbineDomain::isInsideGlobalDomain(const Vector3<real_t> & pt, const real_t& radius) const {
            return isInsideAABB(pt, radius, AABB(blockForest_->getDomain()));
        }


    } // domain

} // turbine_core 

#endif // TURBINECORE_TURBINEDOMAIN_T
