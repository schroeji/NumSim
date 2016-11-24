#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------
class Geometry {
public:
	/// Constructs a default geometry:
        // driven cavity with 128 x 128 grid, no-slip boundary conditions
        // as shown below
        //
        //     u=1, v=0
        //   -------------
        //   |           |
        //u=0|           |u=0
        //v=0|           |v=0
        //   |           |
        //   |           |
        //   -------------
        //     u=0, v=0
	Geometry ();
	Geometry (const Communicator* comm);

	/// Loads a geometry from a file
	void Load (const char* file);

	/// Returns the number of cells in each dimension
	const multi_index_t&	Size	() const;
	/// Returns the total number of cells in each dimension
	const multi_index_t&	TotalSize	() const;
	/// Returns the length of the domain
	const multi_real_t&		Length	() const;
	/// Returns the total length of the domain
	const multi_real_t&		TotalLength	() const;
	/// Returns the meshwidth
	const multi_real_t&		Mesh	() const;

    /// Updates the velocity field u
    void Update_U (Grid *u) const;
    /// Updates the velocity field v
    void Update_V (Grid *v) const;
    /// Updates the pressure field p
    void Update_P (Grid *p) const;
private:
	const Communicator* _comm;

	multi_index_t	_size;
	multi_index_t	_bsize;
	multi_real_t	_length;
	multi_real_t	_blength;
	multi_real_t	_h;

	multi_real_t	_velocity;
	real_t			_pressure;
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
