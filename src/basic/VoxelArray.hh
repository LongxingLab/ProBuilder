#ifndef INCLUDED_basic_VoxelArray_hh
#define INCLUDED_basic_VoxelArray_hh


#include <vector>
#include <array>
#include <iostream>

#include "basic/types.hh"

namespace basic {

using basic::Real;
using basic::Size;

template< class _Float/*=float*/, class _Value/*=float*/ >
class _VoxelArray {
public:
	typedef _VoxelArray<_Float,_Value> THIS;
	static size_t const DIM = 3;
	typedef _Float Float;
	typedef _Value Value;
	typedef std::array< Size, DIM > Indices;
	typedef std::array< Float, DIM > Bounds;

	typedef typename std::vector<Value>::reference data_reference_type;
	typedef typename std::vector<Value>::const_reference const_data_reference_type;
	// typedef util::SimpleArray<DIM,typename BASE::size_type> Indices;
	// typedef util::SimpleArray<DIM,Float> Bounds;

private:
	// All the things that boost::multi_array used to do...

	void
	resize( Indices const & indices ) {

		data_.clear();
		Size the_size = 1;
		for ( Size i = 0; i < DIM; i++ ) {
			the_size *= indices[i];
		}
		data_.resize( the_size );

		shape_ = indices;
	}


public:
	data_reference_type
	operator() ( Indices const & indices ) {
		Size index = 0;
		for ( Size i = 0; i < DIM; ++i ) {
			index *= shape_[i];
			index += indices[i];
		}
		return data_[ index ];
	}

	const_data_reference_type
	operator() ( Indices const & indices ) const {
		Size index = 0;
		for ( Size i = 0; i < DIM; ++i ) {
			index *= shape_[i];
			index += indices[i];
		}
		return data_[ index ];
	}


public:
	std::vector<Value> &
	data( ) {
		return data_;
	}

	std::vector<Value> const &
	data( ) const {
		return data_;
	}

	Size
	num_elements( ) const {
		return data_.size();
	}

	Indices const &
	shape( ) const {
		return shape_;
	}

public:


	_VoxelArray() {}

	template<class F1,class F2, class F3>
	_VoxelArray(F1 const & lb, F2 const & ub, F3 const & cs ) {

		for ( Size i = 0; i < DIM; i++ ) {
			lb_[ i ] = lb[ i ];
			ub_[ i ] = ub[ i ];
			cs_[ i ] = cs[ i ];
		}

		Indices extents = floats_to_index(ub_);
		// std::cout << extents << std::endl;
		for ( Size i = 0; i < DIM; i++ ) extents[i] += 1; // pad by one
		this->resize(extents);
	}

	template<class Floats> Indices floats_to_index(Floats const & f) const {
		Indices ind;
		for ( Size i = 0; i < DIM; ++i ) {
			Float tmp = ((f[i]-lb_[i])/cs_[i]);
			// assert(tmp >= 0.0);
			ind[i] = tmp;
		}
		// std::cout << "floats_to_index " << f << " " << ind << std::endl;
		return ind;
	}

	Bounds indices_to_center( Indices const & idx ) const {
		Bounds c;
		for ( Size i = 0; i < DIM; ++i ) {
			c[i] = (idx[i]+0.5)*cs_[i] + lb_[i];
		}
		return c;
	}

	template<class Floats>
	const_data_reference_type
	operator[](Floats const & floats) const { return this->operator()(floats_to_index(floats)); }

	template<class Floats>
	data_reference_type
	operator[](Floats const & floats){ return this->operator()(floats_to_index(floats)); }

	Value at( Float f, Float g, Float h ) const {
		Indices idx = floats_to_index( Bounds{ { f, g, h } } );
		if ( idx[0] < this->shape()[0] && idx[1] < this->shape()[1] && idx[2] < this->shape()[2] ) {
			return this->operator()(idx);
		} else return Value(0);
	}

	template<class V,typename T=_Value,typename std::enable_if<!std::is_same<T, std::pair<Size,Size>>::value, std::size_t>::type = 0>
	Value at( V const & v ) const {
		Indices idx = floats_to_index( Bounds{ { v[0], v[1], v[2] } } );
		if ( idx[0] < this->shape()[0] && idx[1] < this->shape()[1] && idx[2] < this->shape()[2] && idx[0]>=0 && idx[1]>=0 && idx[2]>=0 ) {
			return _Value(this->operator()(idx));
		} else return _Value(0);
	}
	template<class V,typename T=_Value,typename std::enable_if<std::is_same<T, std::pair<Size,Size>>::value, std::size_t>::type = 0>
	Value at( V const & v ) const {
		Indices idx = floats_to_index( Bounds{ { v[0], v[1], v[2] } } );
		if ( idx[0] < this->shape()[0] && idx[1] < this->shape()[1] && idx[2] < this->shape()[2] && idx[0]>=0 && idx[1]>=0 && idx[2]>=0 ) {
			return _Value(this->operator()(idx));
		} else return _Value(0,0);
	}

	// void write(std::ostream & out) const {
	//  out.write( (char const*)&lb_, sizeof(Bounds) );
	//  out.write( (char const*)&ub_, sizeof(Bounds) );
	//  out.write( (char const*)&cs_, sizeof(Bounds) );
	//  for(size_t i = 0; i < DIM; ++i) out.write( (char const*)&(this->shape()[i]), sizeof() );
	//  out.write( (char const*)this->data(), this->num_elements()*sizeof(Float) );
	// }
	// void read(std::istream & in){
	//  in.read( (char*)&lb_, sizeof(Bounds) );
	//  in.read( (char*)&ub_, sizeof(Bounds) );
	//  in.read( (char*)&cs_, sizeof(Bounds) );
	//  in.read( (char*)this->data(), this->num_elements()*sizeof(Float) );
	// }
	bool operator==(THIS const & o) const {
		return lb_==o.lb_ && ub_==o.ub_ && cs_==o.cs_ && shape_ == o.shape_ && data_ == o.data_;
	}

	void save( std::ostream & out ) const {
		static_assert( std::is_pod<Float>::type::value, "_VoxelArray type is not POD." );
		out.write( (char*)&lb_, sizeof(Bounds) );
		out.write( (char*)&ub_, sizeof(Bounds) );
		out.write( (char*)&cs_, sizeof(Bounds) );
		for ( size_t i = 0; i < DIM; ++i ) {
			out.write( (char*)(&(this->shape()[i])), sizeof(Size) );
		}
		for ( size_t i = 0; i < this->num_elements(); ++i ) out.write( (char*)(&(this->data()[i])), sizeof(Float) );
	}
	void load( std::istream & in ){
		static_assert( std::is_pod<Float>::type::value, "_VoxelArray type is not POD." );
		assert( in.good() );
		in.read( (char*)&lb_, sizeof(Bounds) );
		assert( in.good() );
		in.read( (char*)&ub_, sizeof(Bounds) );
		assert( in.good() );
		in.read( (char*)&cs_, sizeof(Bounds) );
		assert( in.good() );
		// todo: should I check these against the c'tor values? if not, should add default ctor?
		Indices extents;
		for ( size_t i = 0; i < DIM; ++i ) {
			in.read( (char*)(&(extents[i])), sizeof(Size) );
		}
		assert( in.good() );
		this->resize(extents);
		for ( size_t i = 0; i < this->num_elements(); ++i ) in.read( (char*)(&(this->data()[i])), sizeof(Float) );
		assert( in.good() );
	}

	// only for DIM == 3
	// to visualize the pdb grid
	void visualize ( std::ostream & out ) const {
		// only for DIM == 3
		static_assert( std::is_pod<Float>::type::value, "_VoxelArray type is not POD." );

		char buf[128];
		Size anum(1), rnum(1);

		for(Size xx=0; xx<this->shape()[0]; ++xx) {
			for(Size yy=0; yy<this->shape()[1]; ++yy) {
				for(Size zz=0; zz<this->shape()[2]; ++zz) {
					Indices idx = Indices{ { xx, yy, zz } };
					if(0 != this->operator()(idx)) {
						Bounds xyz = this->indices_to_center(idx);
						snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
							"HETATM",
							anum++,
							"BURR",
							"BUR",
							'B',
							rnum++,
							xyz[0],xyz[1],xyz[2],
							1.0,
							1.0,
							"B"
						);

						out << buf;

						rnum %= 10000;
						anum %= 100000;
					}
				}
			}
		}
	}

private:
	Bounds lb_,ub_,cs_;

	// Multi_array stuff
	std::vector<Value> data_;
	Indices shape_;

};

template< class F, class V >
std::ostream & operator << ( std::ostream & out, _VoxelArray<F,V> const & v ){
	out << "_VoxelArray( lb: " << v.lb_ << " ub: " << v.ub_ << " cs: " << v.cs_ << " nelem: " << v.num_elements() << " sizeof_val: " << sizeof(V) << " )";
	return out;
}
// VoxelArray
typedef _VoxelArray<Real,bool> VoxelArray;
typedef _VoxelArray<Real,std::pair<Size,Size>> VoxelIndicateArray;
typedef std::shared_ptr< _VoxelArray<Real,bool> > VoxelArrayOP;
typedef std::shared_ptr< _VoxelArray<Real,std::pair<Size,Size>> > VoxelIndicateArrayOP;


}

#endif
