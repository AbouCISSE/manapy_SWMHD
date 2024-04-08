#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include <stdlib.h>
#include "ndarrays.h"
#include <stdint.h>
#include "cwrapper_ndarrays.h"


void bind_c_explicitscheme_dissipative(int64_t n0_wx_face, double *wx_face, int64_t n0_wy_face, double *wy_face, int64_t n0_wz_face, double *wz_face, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_namef, int64_t *namef, int64_t n0_dissip_w, double *dissip_w, double Dxx, double Dyy, double Dzz);
void bind_c_compute_upwind_flux(double w_l, double w_r, double u_face, double v_face, double w_face, int64_t n0_normal, double *normal, int64_t n0_flux_w, double *flux_w);
void bind_c_explicitscheme_convective_2d(int64_t n0_rez_w, double *rez_w, int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_u_face, double *u_face, int64_t n0_v_face, double *v_face, int64_t n0_w_face, double *w_face, int64_t n0_w_x, double *w_x, int64_t n0_w_y, double *w_y, int64_t n0_w_z, double *w_z, int64_t n0_wx_halo, double *wx_halo, int64_t n0_wy_halo, double *wy_halo, int64_t n0_wz_halo, double *wz_halo, int64_t n0_psi, double *psi, int64_t n0_psi_halo, double *psi_halo, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_centerg, int64_t n1_centerg, double *centerg, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_mesuref, double *mesuref, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_halofid, int64_t *halofid, int64_t n0_name, int64_t *name, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_boundaryfaces, int64_t *boundaryfaces, int64_t n0_periodicboundaryfaces, int64_t *periodicboundaryfaces, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t order);
void bind_c_explicitscheme_convective_3d(int64_t n0_rez_w, double *rez_w, int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_u_face, double *u_face, int64_t n0_v_face, double *v_face, int64_t n0_w_face, double *w_face, int64_t n0_w_x, double *w_x, int64_t n0_w_y, double *w_y, int64_t n0_w_z, double *w_z, int64_t n0_wx_halo, double *wx_halo, int64_t n0_wy_halo, double *wy_halo, int64_t n0_wz_halo, double *wz_halo, int64_t n0_psi, double *psi, int64_t n0_psi_halo, double *psi_halo, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_centerg, int64_t n1_centerg, double *centerg, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_mesuref, double *mesuref, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_halofid, int64_t *halofid, int64_t n0_name, int64_t *name, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_boundaryfaces, int64_t *boundaryfaces, int64_t n0_periodicboundaryfaces, int64_t *periodicboundaryfaces, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t order);

/*........................................*/


/*........................................*/

/*........................................*/
PyObject *explicitscheme_dissipative_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray wx_face = {.shape = NULL};
    t_ndarray wy_face = {.shape = NULL};
    t_ndarray wz_face = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray dissip_w = {.shape = NULL};
    double Dxx;
    double Dyy;
    double Dzz;
    PyArrayObject *wx_face_tmp;
    PyArrayObject *wy_face_tmp;
    PyArrayObject *wz_face_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *dissip_w_tmp;
    PyObject *Dxx_tmp;
    PyObject *Dyy_tmp;
    PyObject *Dzz_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "wx_face",
        "wy_face",
        "wz_face",
        "cellidf",
        "normalf",
        "namef",
        "dissip_w",
        "Dxx",
        "Dyy",
        "Dzz",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!OOO", kwlist, &PyArray_Type, &wx_face_tmp, &PyArray_Type, &wy_face_tmp, &PyArray_Type, &wz_face_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &dissip_w_tmp, &Dxx_tmp, &Dyy_tmp, &Dzz_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(wx_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wx_face = pyarray_to_ndarray(wx_face_tmp);
    }
    if (!pyarray_check(wy_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wy_face = pyarray_to_ndarray(wy_face_tmp);
    }
    if (!pyarray_check(wz_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wz_face = pyarray_to_ndarray(wz_face_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(dissip_w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dissip_w = pyarray_to_ndarray(dissip_w_tmp);
    }
    if (PyIs_NativeFloat(Dxx_tmp))
    {
        Dxx = PyDouble_to_Double(Dxx_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(Dyy_tmp))
    {
        Dyy = PyDouble_to_Double(Dyy_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(Dzz_tmp))
    {
        Dzz = PyDouble_to_Double(Dzz_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    bind_c_explicitscheme_dissipative(nd_ndim(&wx_face, 0), nd_data(&wx_face), nd_ndim(&wy_face, 0), nd_data(&wy_face), nd_ndim(&wz_face, 0), nd_data(&wz_face), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&dissip_w, 0), nd_data(&dissip_w), Dxx, Dyy, Dzz);
    result = Py_BuildValue("");
    free_pointer(wx_face);
    free_pointer(wy_face);
    free_pointer(wz_face);
    free_pointer(cellidf);
    free_pointer(normalf);
    free_pointer(namef);
    free_pointer(dissip_w);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *compute_upwind_flux_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    double w_l;
    double w_r;
    double u_face;
    double v_face;
    double w_face;
    t_ndarray normal = {.shape = NULL};
    t_ndarray flux_w = {.shape = NULL};
    PyObject *w_l_tmp;
    PyObject *w_r_tmp;
    PyObject *u_face_tmp;
    PyObject *v_face_tmp;
    PyObject *w_face_tmp;
    PyArrayObject *normal_tmp;
    PyArrayObject *flux_w_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_l",
        "w_r",
        "u_face",
        "v_face",
        "w_face",
        "normal",
        "flux_w",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOO!O!", kwlist, &w_l_tmp, &w_r_tmp, &u_face_tmp, &v_face_tmp, &w_face_tmp, &PyArray_Type, &normal_tmp, &PyArray_Type, &flux_w_tmp))
    {
        return NULL;
    }
    if (PyIs_NativeFloat(w_l_tmp))
    {
        w_l = PyDouble_to_Double(w_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(w_r_tmp))
    {
        w_r = PyDouble_to_Double(w_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(u_face_tmp))
    {
        u_face = PyDouble_to_Double(u_face_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(v_face_tmp))
    {
        v_face = PyDouble_to_Double(v_face_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(w_face_tmp))
    {
        w_face = PyDouble_to_Double(w_face_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(normal_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        normal = pyarray_to_ndarray(normal_tmp);
    }
    if (!pyarray_check(flux_w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        flux_w = pyarray_to_ndarray(flux_w_tmp);
    }
    bind_c_compute_upwind_flux(w_l, w_r, u_face, v_face, w_face, nd_ndim(&normal, 0), nd_data(&normal), nd_ndim(&flux_w, 0), nd_data(&flux_w));
    result = Py_BuildValue("");
    free_pointer(normal);
    free_pointer(flux_w);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *explicitscheme_convective_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray rez_w = {.shape = NULL};
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray u_face = {.shape = NULL};
    t_ndarray v_face = {.shape = NULL};
    t_ndarray w_face = {.shape = NULL};
    t_ndarray w_x = {.shape = NULL};
    t_ndarray w_y = {.shape = NULL};
    t_ndarray w_z = {.shape = NULL};
    t_ndarray wx_halo = {.shape = NULL};
    t_ndarray wy_halo = {.shape = NULL};
    t_ndarray wz_halo = {.shape = NULL};
    t_ndarray psi = {.shape = NULL};
    t_ndarray psi_halo = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray centerg = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray mesuref = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray name = {.shape = NULL};
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray boundaryfaces = {.shape = NULL};
    t_ndarray periodicboundaryfaces = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t order;
    PyArrayObject *rez_w_tmp;
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *u_face_tmp;
    PyArrayObject *v_face_tmp;
    PyArrayObject *w_face_tmp;
    PyArrayObject *w_x_tmp;
    PyArrayObject *w_y_tmp;
    PyArrayObject *w_z_tmp;
    PyArrayObject *wx_halo_tmp;
    PyArrayObject *wy_halo_tmp;
    PyArrayObject *wz_halo_tmp;
    PyArrayObject *psi_tmp;
    PyArrayObject *psi_halo_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerf_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *centerg_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *mesuref_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *name_tmp;
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *boundaryfaces_tmp;
    PyArrayObject *periodicboundaryfaces_tmp;
    PyArrayObject *shift_tmp;
    PyObject *order_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "rez_w",
        "w_c",
        "w_ghost",
        "w_halo",
        "u_face",
        "v_face",
        "w_face",
        "w_x",
        "w_y",
        "w_z",
        "wx_halo",
        "wy_halo",
        "wz_halo",
        "psi",
        "psi_halo",
        "centerc",
        "centerf",
        "centerh",
        "centerg",
        "cellidf",
        "mesuref",
        "normalf",
        "halofid",
        "name",
        "innerfaces",
        "halofaces",
        "boundaryfaces",
        "periodicboundaryfaces",
        "shift",
        "order",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O", kwlist, &PyArray_Type, &rez_w_tmp, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &u_face_tmp, &PyArray_Type, &v_face_tmp, &PyArray_Type, &w_face_tmp, &PyArray_Type, &w_x_tmp, &PyArray_Type, &w_y_tmp, &PyArray_Type, &w_z_tmp, &PyArray_Type, &wx_halo_tmp, &PyArray_Type, &wy_halo_tmp, &PyArray_Type, &wz_halo_tmp, &PyArray_Type, &psi_tmp, &PyArray_Type, &psi_halo_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerf_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &centerg_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &mesuref_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &name_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &boundaryfaces_tmp, &PyArray_Type, &periodicboundaryfaces_tmp, &PyArray_Type, &shift_tmp, &order_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(rez_w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_w = pyarray_to_ndarray(rez_w_tmp);
    }
    if (!pyarray_check(w_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_c = pyarray_to_ndarray(w_c_tmp);
    }
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(u_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u_face = pyarray_to_ndarray(u_face_tmp);
    }
    if (!pyarray_check(v_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_face = pyarray_to_ndarray(v_face_tmp);
    }
    if (!pyarray_check(w_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_face = pyarray_to_ndarray(w_face_tmp);
    }
    if (!pyarray_check(w_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_x = pyarray_to_ndarray(w_x_tmp);
    }
    if (!pyarray_check(w_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_y = pyarray_to_ndarray(w_y_tmp);
    }
    if (!pyarray_check(w_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_z = pyarray_to_ndarray(w_z_tmp);
    }
    if (!pyarray_check(wx_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wx_halo = pyarray_to_ndarray(wx_halo_tmp);
    }
    if (!pyarray_check(wy_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wy_halo = pyarray_to_ndarray(wy_halo_tmp);
    }
    if (!pyarray_check(wz_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wz_halo = pyarray_to_ndarray(wz_halo_tmp);
    }
    if (!pyarray_check(psi_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi = pyarray_to_ndarray(psi_tmp);
    }
    if (!pyarray_check(psi_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi_halo = pyarray_to_ndarray(psi_halo_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = pyarray_to_ndarray(centerf_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(centerg_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerg = pyarray_to_ndarray(centerg_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(mesuref_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesuref = pyarray_to_ndarray(mesuref_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(name_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        name = pyarray_to_ndarray(name_tmp);
    }
    if (!pyarray_check(innerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        innerfaces = pyarray_to_ndarray(innerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(boundaryfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        boundaryfaces = pyarray_to_ndarray(boundaryfaces_tmp);
    }
    if (!pyarray_check(periodicboundaryfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        periodicboundaryfaces = pyarray_to_ndarray(periodicboundaryfaces_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(order_tmp))
    {
        order = PyInt64_to_Int64(order_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_explicitscheme_convective_2d(nd_ndim(&rez_w, 0), nd_data(&rez_w), nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&u_face, 0), nd_data(&u_face), nd_ndim(&v_face, 0), nd_data(&v_face), nd_ndim(&w_face, 0), nd_data(&w_face), nd_ndim(&w_x, 0), nd_data(&w_x), nd_ndim(&w_y, 0), nd_data(&w_y), nd_ndim(&w_z, 0), nd_data(&w_z), nd_ndim(&wx_halo, 0), nd_data(&wx_halo), nd_ndim(&wy_halo, 0), nd_data(&wy_halo), nd_ndim(&wz_halo, 0), nd_data(&wz_halo), nd_ndim(&psi, 0), nd_data(&psi), nd_ndim(&psi_halo, 0), nd_data(&psi_halo), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&centerg, 0), nd_ndim(&centerg, 1), nd_data(&centerg), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&mesuref, 0), nd_data(&mesuref), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&name, 0), nd_data(&name), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&boundaryfaces, 0), nd_data(&boundaryfaces), nd_ndim(&periodicboundaryfaces, 0), nd_data(&periodicboundaryfaces), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), order);
    result = Py_BuildValue("");
    free_pointer(rez_w);
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(u_face);
    free_pointer(v_face);
    free_pointer(w_face);
    free_pointer(w_x);
    free_pointer(w_y);
    free_pointer(w_z);
    free_pointer(wx_halo);
    free_pointer(wy_halo);
    free_pointer(wz_halo);
    free_pointer(psi);
    free_pointer(psi_halo);
    free_pointer(centerc);
    free_pointer(centerf);
    free_pointer(centerh);
    free_pointer(centerg);
    free_pointer(cellidf);
    free_pointer(mesuref);
    free_pointer(normalf);
    free_pointer(halofid);
    free_pointer(name);
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(boundaryfaces);
    free_pointer(periodicboundaryfaces);
    free_pointer(shift);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *explicitscheme_convective_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray rez_w = {.shape = NULL};
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray u_face = {.shape = NULL};
    t_ndarray v_face = {.shape = NULL};
    t_ndarray w_face = {.shape = NULL};
    t_ndarray w_x = {.shape = NULL};
    t_ndarray w_y = {.shape = NULL};
    t_ndarray w_z = {.shape = NULL};
    t_ndarray wx_halo = {.shape = NULL};
    t_ndarray wy_halo = {.shape = NULL};
    t_ndarray wz_halo = {.shape = NULL};
    t_ndarray psi = {.shape = NULL};
    t_ndarray psi_halo = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray centerg = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray mesuref = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray name = {.shape = NULL};
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray boundaryfaces = {.shape = NULL};
    t_ndarray periodicboundaryfaces = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t order;
    PyArrayObject *rez_w_tmp;
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *u_face_tmp;
    PyArrayObject *v_face_tmp;
    PyArrayObject *w_face_tmp;
    PyArrayObject *w_x_tmp;
    PyArrayObject *w_y_tmp;
    PyArrayObject *w_z_tmp;
    PyArrayObject *wx_halo_tmp;
    PyArrayObject *wy_halo_tmp;
    PyArrayObject *wz_halo_tmp;
    PyArrayObject *psi_tmp;
    PyArrayObject *psi_halo_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerf_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *centerg_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *mesuref_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *name_tmp;
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *boundaryfaces_tmp;
    PyArrayObject *periodicboundaryfaces_tmp;
    PyArrayObject *shift_tmp;
    PyObject *order_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "rez_w",
        "w_c",
        "w_ghost",
        "w_halo",
        "u_face",
        "v_face",
        "w_face",
        "w_x",
        "w_y",
        "w_z",
        "wx_halo",
        "wy_halo",
        "wz_halo",
        "psi",
        "psi_halo",
        "centerc",
        "centerf",
        "centerh",
        "centerg",
        "cellidf",
        "mesuref",
        "normalf",
        "halofid",
        "name",
        "innerfaces",
        "halofaces",
        "boundaryfaces",
        "periodicboundaryfaces",
        "shift",
        "order",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O", kwlist, &PyArray_Type, &rez_w_tmp, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &u_face_tmp, &PyArray_Type, &v_face_tmp, &PyArray_Type, &w_face_tmp, &PyArray_Type, &w_x_tmp, &PyArray_Type, &w_y_tmp, &PyArray_Type, &w_z_tmp, &PyArray_Type, &wx_halo_tmp, &PyArray_Type, &wy_halo_tmp, &PyArray_Type, &wz_halo_tmp, &PyArray_Type, &psi_tmp, &PyArray_Type, &psi_halo_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerf_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &centerg_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &mesuref_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &name_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &boundaryfaces_tmp, &PyArray_Type, &periodicboundaryfaces_tmp, &PyArray_Type, &shift_tmp, &order_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(rez_w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_w = pyarray_to_ndarray(rez_w_tmp);
    }
    if (!pyarray_check(w_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_c = pyarray_to_ndarray(w_c_tmp);
    }
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(u_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u_face = pyarray_to_ndarray(u_face_tmp);
    }
    if (!pyarray_check(v_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_face = pyarray_to_ndarray(v_face_tmp);
    }
    if (!pyarray_check(w_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_face = pyarray_to_ndarray(w_face_tmp);
    }
    if (!pyarray_check(w_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_x = pyarray_to_ndarray(w_x_tmp);
    }
    if (!pyarray_check(w_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_y = pyarray_to_ndarray(w_y_tmp);
    }
    if (!pyarray_check(w_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_z = pyarray_to_ndarray(w_z_tmp);
    }
    if (!pyarray_check(wx_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wx_halo = pyarray_to_ndarray(wx_halo_tmp);
    }
    if (!pyarray_check(wy_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wy_halo = pyarray_to_ndarray(wy_halo_tmp);
    }
    if (!pyarray_check(wz_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wz_halo = pyarray_to_ndarray(wz_halo_tmp);
    }
    if (!pyarray_check(psi_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi = pyarray_to_ndarray(psi_tmp);
    }
    if (!pyarray_check(psi_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi_halo = pyarray_to_ndarray(psi_halo_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = pyarray_to_ndarray(centerf_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(centerg_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerg = pyarray_to_ndarray(centerg_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(mesuref_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesuref = pyarray_to_ndarray(mesuref_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(name_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        name = pyarray_to_ndarray(name_tmp);
    }
    if (!pyarray_check(innerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        innerfaces = pyarray_to_ndarray(innerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(boundaryfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        boundaryfaces = pyarray_to_ndarray(boundaryfaces_tmp);
    }
    if (!pyarray_check(periodicboundaryfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        periodicboundaryfaces = pyarray_to_ndarray(periodicboundaryfaces_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(order_tmp))
    {
        order = PyInt64_to_Int64(order_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_explicitscheme_convective_3d(nd_ndim(&rez_w, 0), nd_data(&rez_w), nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&u_face, 0), nd_data(&u_face), nd_ndim(&v_face, 0), nd_data(&v_face), nd_ndim(&w_face, 0), nd_data(&w_face), nd_ndim(&w_x, 0), nd_data(&w_x), nd_ndim(&w_y, 0), nd_data(&w_y), nd_ndim(&w_z, 0), nd_data(&w_z), nd_ndim(&wx_halo, 0), nd_data(&wx_halo), nd_ndim(&wy_halo, 0), nd_data(&wy_halo), nd_ndim(&wz_halo, 0), nd_data(&wz_halo), nd_ndim(&psi, 0), nd_data(&psi), nd_ndim(&psi_halo, 0), nd_data(&psi_halo), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&centerg, 0), nd_ndim(&centerg, 1), nd_data(&centerg), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&mesuref, 0), nd_data(&mesuref), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&name, 0), nd_data(&name), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&boundaryfaces, 0), nd_data(&boundaryfaces), nd_ndim(&periodicboundaryfaces, 0), nd_data(&periodicboundaryfaces), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), order);
    result = Py_BuildValue("");
    free_pointer(rez_w);
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(u_face);
    free_pointer(v_face);
    free_pointer(w_face);
    free_pointer(w_x);
    free_pointer(w_y);
    free_pointer(w_z);
    free_pointer(wx_halo);
    free_pointer(wy_halo);
    free_pointer(wz_halo);
    free_pointer(psi);
    free_pointer(psi_halo);
    free_pointer(centerc);
    free_pointer(centerf);
    free_pointer(centerh);
    free_pointer(centerg);
    free_pointer(cellidf);
    free_pointer(mesuref);
    free_pointer(normalf);
    free_pointer(halofid);
    free_pointer(name);
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(boundaryfaces);
    free_pointer(periodicboundaryfaces);
    free_pointer(shift);
    return result;
}
/*........................................*/

static int exec_func(PyObject* m)
{
    return 0;
}

/*........................................*/

static PyMethodDef pyccel_fvm_methods[] = {
    {
        "explicitscheme_dissipative",
        (PyCFunction)explicitscheme_dissipative_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "compute_upwind_flux",
        (PyCFunction)compute_upwind_flux_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "explicitscheme_convective_2d",
        (PyCFunction)explicitscheme_convective_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "explicitscheme_convective_3d",
        (PyCFunction)explicitscheme_convective_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static PyModuleDef_Slot pyccel_fvm_slots[] = {
    {Py_mod_exec, exec_func},
    {0, NULL},
};

/*........................................*/

static struct PyModuleDef pyccel_fvm_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "pyccel_fvm",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    pyccel_fvm_methods,
    pyccel_fvm_slots
};

/*........................................*/

PyMODINIT_FUNC PyInit_pyccel_fvm(void)
{
    import_array();
    return PyModuleDef_Init(&pyccel_fvm_module);
}
