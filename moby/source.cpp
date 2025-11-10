#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <vector>
#include <cmath>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Constant.h>


namespace py = pybind11;

class ConeBeam: public dolfin::Expression {
    public:

    ConeBeam(int dim): dolfin::Expression(dim), _dim(dim), _xs(dim), _d(dim) {}
    
    
    void set(py::array_t<double>& x, py::array_t<double>& d, double theta_max, double intensity, double mua)
    {
        auto xx = x.unchecked<1>(); 
        auto dd = d.unchecked<1>();
        double norm_d = 0.;
        for(int i(0); i < _dim; ++i)
        {
            _xs[i] = xx[i];
            _d[i]  = dd[i];
            norm_d += _d[i]*_d[i];
        }
        norm_d = std::sqrt(norm_d);
        for(int i(0); i < _dim; ++i)
            _d[i] /= norm_d;

        _theta_max = theta_max;
        _cos_theta_max = std::cos(theta_max);
        _intensity = intensity;
        _mua = mua;
    }

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> point, const ufc::cell& cell) const final {
        values = point - _xs;
        double distance2 = values.dot(values);
        double distance = std::sqrt(distance2);
        double cos_theta = values.dot(_d)/distance;
        if (cos_theta < _cos_theta_max)
        values *= 0.;
        else
        values *= _intensity/(4.* DOLFIN_PI*distance2*distance)*std::exp(-_mua*distance)*cos_theta;
    }

    private:
    int _dim;
    Eigen::VectorXd _xs;
    Eigen::VectorXd _d;
    double _theta_max;
    double _cos_theta_max;
    double _intensity;
    double _mua;
};


class SlitBeam: public dolfin::Expression {
    public:

    SlitBeam(int dim): dolfin::Expression(dim), _dim(dim), _xs(dim), _d(dim) {}
    
    
    void set(py::array_t<double>& x, py::array_t<double>& d, double height, double theta_max, double intensity, double mua)
    {
        auto xx = x.unchecked<1>(); 
        auto dd = d.unchecked<1>();
        for(int i(0); i < _dim; ++i)
        {
            _xs[i] = xx[i];
            _d[i]  = dd[i];
        }

        assert( std::abs(_d[_dim-1]) < 1e-16 ); //The direction of propagation must be in the x-y plane
        assert ( std::abs( _d.dot(_d) - 1. ) < 1e-8 ); //The norm of _d must be 1

        _half_height    = .5*height;
        _theta_max = theta_max;
        _cos_theta_max = std::cos(theta_max);
        _intensity = intensity;
        _mua = mua;
    }

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> point, const ufc::cell& cell) const final {
        values = point - _xs;
        double z_diff = values[_dim-1]; 
        values[_dim-1] = 0.; // Rays only propagate in the x-y plane
        double distance2 = values.dot(values);
        double distance = std::sqrt(distance2);
        double cos_theta_num = values.dot(_d);
        double cos_theta = values.dot(_d)/distance;
        if (( std::abs(z_diff) <=  _half_height ) && (cos_theta >= _cos_theta_max)) //Illumination is non-zero only inside the prism
            // values *= _intensity/(4.* DOLFIN_PI*distance2*distance)*std::exp(-_mua*distance)*cos_theta;
            values *= _intensity/(2.* DOLFIN_PI*distance2)*std::exp(-_mua*distance)*cos_theta;
        else
            values *= 0.;
    }

    private:
    int _dim;            //Geometric dimension (should be 3)
    Eigen::VectorXd _xs; // Center of the slit
    Eigen::VectorXd _d;  // Direction of the slit (must be perpendicular to z-axis)
    double _half_height; // Half height of the slit
    double _theta_max;   // Maximum angle in the x-y plane
    double _cos_theta_max;
    double _intensity;   //Maximum illumination
    double _mua;         // Absorption coefficient units 1/L where L is the units used for the mesh coordinates.
};

class BoundaryNormal : public dolfin::Expression
{
  public:
    BoundaryNormal(std::shared_ptr<dolfin::Mesh> mesh): dolfin::Expression( mesh->geometry().dim() ),
    m_gdim( mesh->geometry().dim() ), m_mesh(mesh){}

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> point,
                                            const ufc::cell &cell) const final
    {

      dolfin::Cell c(*m_mesh,cell.index);
      
      dolfin::MeshEntityIterator fac_it(c, m_gdim-1);
      dolfin::Point normal;

      for(; !fac_it.end(); ++fac_it)
      {
        dolfin::Facet f(*m_mesh,fac_it->index());
        if(f.exterior())
          normal = f.normal();
      }

      values[0] = normal[0];
      values[1] = normal[1];
      if(m_gdim == 3) values[2] = normal[2];
    }

    int gdim() const { return m_gdim; }

  private:
    int m_gdim;
    std::shared_ptr<dolfin::Mesh> m_mesh;
};

class BoundaryQ0Source: public dolfin::Expression
{
    public:
    BoundaryQ0Source(std::shared_ptr<dolfin::Mesh> mesh):
            dolfin::Expression(),
            m_normal(  BoundaryNormal(mesh) ),
            m_S(0)
            {}

    void append(std::shared_ptr<dolfin::Expression> & Si)
    {
        m_S.push_back(Si);
    }

    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> point,
              const ufc::cell &cell) const final
    {
        double out(0.);
        double outi(0.);
        Eigen::VectorXd S(m_normal.gdim());
        Eigen::VectorXd n(m_normal.gdim());
        m_normal.eval(n, point, cell);
        for( auto Si(m_S.begin()); Si != m_S.end(); ++Si)
        {
            (*Si)->eval(S, point, cell);
            outi = -n.dot(S);
            if(outi > 0.) out += outi;
        }

        values[0] = out;
    }

    private:
    BoundaryNormal m_normal;
    std::vector< std::shared_ptr<dolfin::Expression> > m_S;
};


PYBIND11_MODULE(SIGNATURE, m)
    {
    py::class_<ConeBeam, std::shared_ptr<ConeBeam>, dolfin::Expression>
    (m, "ConeBeam", "A cone beam light source in 2 or 3 dimension")
    .def(py::init<int>())
    .def("set", &ConeBeam::set, "x: Cone coordinates; d: Direction of the cone; theta_max: cone angle; intensity: Intensity of beam; mua: absorption coefficient of coupling medium");

    py::class_<SlitBeam, std::shared_ptr<SlitBeam>, dolfin::Expression>
    (m, "SlitBeam", "A slit source in 3 dimension. The slit is assumed to be parallel to the z-axis; photons are emitted perpendicularly to the slit")
    .def(py::init<int>())
    .def("set", &SlitBeam::set, "x: Center of the slit; d: Direction of the cone in the x-y plane; height: Height of slit in z-direction; theta_max: cone angle in the x-y plane; intensity: Intensity of beam; mua: absorption coefficient of coupling medium");

    py::class_<BoundaryNormal, std::shared_ptr<BoundaryNormal>, dolfin::Expression>
    (m, "BoundaryNormal", "Computes the normal vector to boundary cells")
    .def(py::init< std::shared_ptr<dolfin::Mesh> >());

    py::class_<BoundaryQ0Source,  std::shared_ptr<BoundaryQ0Source>, dolfin::Expression>
    (m, "BoundaryQ0Source")
    .def( py::init< std::shared_ptr<dolfin::Mesh> >() )
    .def( "append", &BoundaryQ0Source::append );
    }

