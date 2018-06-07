#ifndef VEC3_HPP
#define VEC3_HPP

#include <iostream>
#include <math.h>

/*
 * Template for 3d vectors
 */

template <class T> class Vec3
{
  private:
    T x, y, z;

  public:
    Vec3()
    {
        x = y = z = 0;
    };

    Vec3(T xValue, T yValue, T zValue)
    {
        x = xValue;
        y = yValue;
        z = zValue;
    }

    void set(const T &xValue, const T &yValue, const T &zValue)
    {
        x = xValue;
        y = yValue;
        z = zValue;
    }

    T getX() const
    {
        return x;
    }
    T getY() const
    {
        return y;
    }
    T getZ() const
    {
        return z;
    }

    void setX(const T &xValue)
    {
        x = xValue;
    }
    void setY(const T &yValue)
    {
        y = yValue;
    }
    void setZ(const T &zValue)
    {
        z = zValue;
    }

    void zero()
    {
        x = y = z = 0;
    }

    void normalise()
    {
        T magnitude = sqrt((x * x) + (y * y) + (z * z));

        if (magnitude != 0)
        {
            x /= magnitude;
            y /= magnitude;
            z /= magnitude;
        }
    }

    static T dotProduct(const Vec3 &vec1, const Vec3 &vec2)
    {
        return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    }

    T dotProduct(const Vec3 &vec) const
    {
        return x * vec.x + y * vec.y + z * vec.z;
    }

    static T getDistance(const Vec3 &v1, const Vec3 &v2)
    {
        T dx = v2.x - v1.x;
        T dy = v2.y - v1.y;
        T dz = v2.z - v1.z;

        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    T abs() const
    {
        return sqrt((x * x) + (y * y) + (z * z));
    }

    T abs2() const
    {
        return (x * x) + (y * y) + (z * z);
    }

    Vec3 operator+(const Vec3 &vector) const
    {
        return Vec3<T>(x + vector.x, y + vector.y, z + vector.z);
    }

    void operator+=(const Vec3 &vector)
    {
        x += vector.x;
        y += vector.y;
        z += vector.z;
    }

    Vec3 operator-(const Vec3 &vector) const
    {
        return Vec3<T>(x - vector.x, y - vector.y, z - vector.z);
    }

    void operator-=(const Vec3 &vector)
    {
        x -= vector.x;
        y -= vector.y;
        z -= vector.z;
    }

    Vec3 operator*(const Vec3 &vector) const
    {
        return Vec3<T>(x * vector.x, y * vector.y, z * vector.z);
    }

    Vec3 operator*(const T &value) const
    {
        return Vec3<T>(x * value, y * value, z * value);
    }

    void operator*=(const T &value)
    {
        x *= value;
        y *= value;
        z *= value;
    }

    Vec3 operator/(const T &value) const
    {
        return Vec3<T>(x / value, y / value, z / value);
    }

    void operator/=(const T &value)
    {
        x /= value;
        y /= value;
        z /= value;
    }
};

#endif
