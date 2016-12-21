//
// Some extra functions missing in STL and Boost.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_UTIL_
#define _FLUKE_UTIL_

namespace fluke {

    /// A variant of the STL function \c std::for_each. 
    /// \pre \c first2 is not in the range [\c first, \c last)
    template< class In, class In2, class BinOp > BinOp
    for_each( In first, In last, In2 first2, BinOp op ) {
        while( first != last ) op( *first++, *first2++ );
        return op;
    }
    
    /// Copy if a predicate holds. A function clearly missing from the STL.
    /// \pre \c res is not in the range [\c first, \c last)
    template< class In, class Out, class Pred > Out
    copy_if( In first, In last, Out res, Pred p ) {
        while( first != last ) {
            if( p( *first ) ) *res++ = *first;
            ++first;
        }
        return res;
    }

    /// Indirect erase. Normally the STL function family \c erase deletes 
    /// only the object in the container. In \b fluke most of the containers 
    /// hold pointers to objects. \c smart_erase erases the object pointed to
    /// as well.
    template< class Container, class For > For
    smart_erase( Container &x, For first ) {
        delete *first;
        return x.erase( first );
    }

    /// Indirect range erase. Both the range and the objects pointed to are
    /// erased.
    /// See also smart_erase( Container, For )
    template< class Container, class For > For
    smart_erase( Container &x, For first, For last ) {
        For i = first;
        while( i != last ) {
            delete *i;
            ++i;
        }
        return x.erase( first, last );
    }
    
    /// Numerically safe floating-point comparison.
    template< class T > bool
    close_to( T f, T g ) {
        return std::abs( f - g ) < std::numeric_limits< T >::epsilon();
    }

    /// Mean
    template< class In > double
    mean( In first, In last ) {
        double mm = 0.0;
        mm = std::accumulate( first, last, mm );
        return mm / std::distance( first, last );
    }

    /// Variance
    template< class In > double
    variance( In first, In last, double mean ) {
        double aux = std::distance( first, last );
        double bux = 0.0;
        while( first != last ) {
            bux += ( *first - mean ) * ( *first - mean );
            ++first;
        }
        return bux / aux;
    }

    /// Retrieve random element
    template< class For, class RandomGenerator > For
    random_element( For first, For last, RandomGenerator &rangen ) {
        int aux = std::distance( first, last );
        return boost::next( first, rangen( aux ) );
    }
}
#endif

