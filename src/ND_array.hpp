#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <cstddef>
#include <iostream>
#include <vector>
#include <omp.h>


#define PAR_SIZE 1024



template <size_t firstDim, size_t... RestDims>
class Array
{
public:

    static constexpr size_t N = sizeof...(RestDims) + 1;
    static constexpr size_t length = firstDim*(RestDims * ...);
    static constexpr size_t Dims[] = {firstDim, RestDims...};
    
    double* data_;

    bool is_original;


    // Base constructor
    Array()
    {
        data_ = new double[length];
        
        std::fill_n(data_,length, 0.0);

        is_original =  true;
    }

    // copy constructor
    Array(const Array<firstDim, RestDims...>& other)
    {
        data_ = other.data_;
        is_original =  false;
    }
    
    // Constructor from pointer
    Array(double* p, bool is_or = false)
    {
        data_ = p;
        is_original = is_or;
    }



    const Array<firstDim, RestDims...>& operator=(const Array<firstDim, RestDims...>& other)
    {
        //data_ = other.data_;
        std::copy(other.data_, other.data_ + length, data_);

        return *this;
    }
    constexpr const Array<firstDim, RestDims...>& operator=(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    
    // constructor from N-1 dimensional array
    constexpr Array(Array<RestDims...> const& slice)
    {
        data_ = new double[length];
        for(size_t i = 0; i<firstDim; i++)
        {
            std::copy(slice.data_, slice.data_ + slice.length, data_+ i*(RestDims*...));
        }
    }

    ~Array() {  if(is_original)    delete[] data_;   }



    // access element
    template <typename... Indices>
    constexpr double& operator()(Indices... indices)
    {
        size_t offset = 0;
        size_t temp[] = {static_cast<size_t>(indices)...};

        for (size_t i = 0; i < N; i++) {
            offset = offset * Dims[i] + temp[i];
        }
        return data_[offset];
    }
    template <typename... Indices>
    constexpr double operator()(Indices... indices) const
    {
        size_t offset = 0;
        size_t temp[] = {static_cast<size_t>(indices)...};

        for (size_t i = 0; i < N; i++) {
            offset = offset * Dims[i] + temp[i];
        }
        return data_[offset];
    }


    // access element
    constexpr Array<RestDims...> operator[](size_t index)
    {
        return Array<RestDims...>(data_+ index * (RestDims * ...));
    }
    constexpr const Array<RestDims...> operator[](size_t index) const
    {
        return Array<RestDims...>(data_+ index * (RestDims * ...));
    }





    constexpr const Array<firstDim, RestDims...>& fill(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    constexpr const Array<firstDim, RestDims...>& fill(const Array<firstDim, RestDims...>& other)
    {
        if(data_ != other.data_)    std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }


    // explicit copy
    constexpr Array<firstDim, RestDims...> copy() const
    {
        Array<firstDim, RestDims...> result;

        std::copy(data_, data_+length, result.data_);

        return result;
    }



    constexpr size_t size(const size_t index = 0) const
    {
        return Dims[index];
    }

    
    // abs
    Array<firstDim, RestDims...> const& abs()
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        return *this;
    }


    // min
    double min() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        return res;
    }
    // max
    double max() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        return res;
    }





    // Arithmetic operations


    // += operator
    const Array<firstDim, RestDims...>& operator+=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator+=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        return *this;
    }

    // -= operator
    const Array<firstDim, RestDims...>& operator-=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator-=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        return *this;
    }

    // *= operator
    const Array<firstDim, RestDims...>& operator*=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator*=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        return *this;
    }

    // /= operator
    const Array<firstDim, RestDims...>& operator/=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator/=(double scalar)
    {
        double inv_scal = 1.0/scalar;

        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        return *this;
    }
};




// ostream
template <size_t firstDim, size_t... RestDims>
std::ostream& operator<<(std::ostream& output, const Array<firstDim, RestDims...>& other)
{
    if constexpr (sizeof...(RestDims) > 0)
    {
        for(size_t i = 0; i<firstDim; i++)
            output<<other[i]<<std::endl;
    }
    else
    {
        for(size_t i = 0; i<firstDim; i++)
            output<<other[i]<<"\t";
    }

    return output;
}





// abs
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> abs(Array<firstDim, RestDims...> M)
{
    Array<firstDim, RestDims...> result = M.copy();

    return result.abs();
}


// min
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> min(Array<firstDim, RestDims...> M)
{
    return M.min();
}
// max
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> max(Array<firstDim, RestDims...> M)
{
    return M.max();
}



// addition
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator+(const Array<firstDim, RestDims...>& M1, const Array<firstDim, RestDims...>& M2)
{
    Array<firstDim, RestDims...> result = M1.copy();
    result+=M2;
    return result;
}
// scalar
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator+(double scalar, const Array<firstDim, RestDims...>& other)
{
    Array<firstDim, RestDims...> result = other.copy();
    result+=scalar;
    return result;
}
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator+(const Array<firstDim, RestDims...>& other, double scalar)
{
    Array<firstDim, RestDims...> result = other.copy();
    result+=scalar;
    return result;
}


// substraction
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator-(const Array<firstDim, RestDims...>& M1, const Array<firstDim, RestDims...>& M2)
{
    Array<firstDim, RestDims...> result = M1.copy();
    result-=M2;
    return result;
}
// scalar
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator-(double scalar, const Array<firstDim, RestDims...>& other)
{
    Array<firstDim, RestDims...> result = other.copy();
    result-=scalar;
    return result;
}
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator-(const Array<firstDim, RestDims...>& other, double scalar)
{
    Array<firstDim, RestDims...> result = other.copy();
    result-=scalar;
    return result;
}


// multiplication
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator*(const Array<firstDim, RestDims...>& M1, const Array<firstDim, RestDims...>& M2)
{
    Array<firstDim, RestDims...> result = M1.copy();
    result*=M2;
    return result;
}
// scalar
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator*(double scalar, const Array<firstDim, RestDims...>& other)
{
    Array<firstDim, RestDims...> result = other.copy();
    result*=scalar;
    return result;
}
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator*(const Array<firstDim, RestDims...>& other, double scalar)
{
    Array<firstDim, RestDims...> result = other.copy();
    result*=scalar;
    return result;
}


// division
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator/(const Array<firstDim, RestDims...>& M1, const Array<firstDim, RestDims...>& M2)
{
    Array<firstDim, RestDims...> result = M1.copy();
    result/=M2;
    return result;
}
// scalar
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator/(double scalar, const Array<firstDim, RestDims...>& other)
{
    Array<firstDim, RestDims...> result;

    result = scalar;

    result/=other;

    return result;
}
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator/(const Array<firstDim, RestDims...>& other, double scalar)
{
    Array<firstDim, RestDims...> result = other.copy();
    result/=scalar;
    return result;
}

// opposite
template <size_t firstDim, size_t... RestDims>
constexpr Array<firstDim, RestDims...> operator-(const Array<firstDim, RestDims...>& other)
{
    Array<firstDim, RestDims...> result;

    result-=other;

    return result;
}




template <size_t Dim>
class Array<Dim>
{
public:

    static constexpr size_t N = 1;
    static constexpr size_t length = Dim;

    double* data_;

    bool is_original;

    
    // Base constructor
    Array()
    {
        data_ = new double[length];
        std::fill_n(data_,length, 0.0);

        is_original =  true;
    }

    // copy constructor
    Array(const Array<Dim> & other)
    {
        data_ = other.data_;
        is_original =  false;
    }
    
    // Constructor from pointer
    Array(double* p, bool is_or = false)
    {
        data_ = p;
        is_original = is_or;
    }
    


    const Array<Dim>& operator=(const Array<Dim>& other)
    {
        //data_ = other.data_;
        std::copy(other.data_, other.data_ + length, data_);

        return *this;
    }
    constexpr const Array<Dim>& operator=(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    
    ~Array() {  if(is_original)    delete[] data_;   }


    // access element
    constexpr double& operator()(size_t index)              { return data_[index]; }
    constexpr double operator()(size_t index) const   { return data_[index]; }
    constexpr double& operator[](size_t index)              { return data_[index]; }
    constexpr double operator[](size_t index) const   { return data_[index]; }





    constexpr const Array<Dim>& fill(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    constexpr const Array<Dim>& fill(const Array<Dim> & other)
    {
        if(data_ != other.data_)    std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }


    // explicit copy
    constexpr Array<Dim> copy() const
    {
        Array<Dim> result;

        std::copy(data_, data_+length, result.data_);

        return result;
    }




    // abs
    Array<Dim> const& abs()
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        return *this;
    }


    // min
    double min() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        return res;
    }
    // max
    double max() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        return res;
    }





    // Arithmetic operations


    // += operator
    const Array<Dim>& operator+=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator+=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        return *this;
    }

    // -= operator
    const Array<Dim>& operator-=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator-=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        return *this;
    }

    // *= operator
    const Array<Dim>& operator*=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator*=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        return *this;
    }

    // /= operator
    const Array<Dim>& operator/=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator/=(double scalar)
    {
        double inv_scal = 1.0/scalar;

        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        return *this;
    }
};



#endif
