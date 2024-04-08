#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include <stdlib.h>
#include "ndarrays.h"
#include "cwrapper_ndarrays.h"
#include <stdint.h>


int bind_c_total_energy(int64_t n0_h_c, double *h_c, int64_t n0_hu_c, double *hu_c, int64_t n0_hv_c, double *hv_c, int64_t n0_hB1_c, double *hB1_c, int64_t n0_hB2_c, double *hB2_c, int64_t n0_Z_c, double *Z_c, double grav, int64_t n0_volumec, double *volumec, double *numt, double *numc, double *numm, double *nump);
void bind_c_explicitscheme_stabilizater(int64_t n0_wx_face, double *wx_face, int64_t n0_wy_face, double *wy_face, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_namef, int64_t *namef, double vepsilon, int64_t n0_dissip_w, double *dissip_w, int64_t n0_h_c, double *h_c);
double bind_c_time_step_swmhd(int64_t n0_h_c, double *h_c, int64_t n0_hu_c, double *hu_c, int64_t n0_hv_c, double *hv_c, int64_t n0_hB1_c, double *hB1_c, int64_t n0_hB2_c, double *hB2_c, double cfl, int64_t n0_normal, int64_t n1_normal, double *normal, int64_t n0_mesure, double *mesure, int64_t n0_volume, double *volume, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid);
void bind_c_update_swmhd(int64_t n0_h_c, double *h_c, int64_t n0_hu_c, double *hu_c, int64_t n0_hv_c, double *hv_c, int64_t n0_hB1_c, double *hB1_c, int64_t n0_hB2_c, double *hB2_c, int64_t n0_PSI_c, double *PSI_c, int64_t n0_Z_c, double *Z_c, int64_t n0_rez_h, double *rez_h, int64_t n0_rez_hu, double *rez_hu, int64_t n0_rez_hv, double *rez_hv, int64_t n0_rez_hB1, double *rez_hB1, int64_t n0_rez_hB2, double *rez_hB2, int64_t n0_rez_PSI, double *rez_PSI, int64_t n0_rez_Z, double *rez_Z, int64_t n0_src_h, double *src_h, int64_t n0_src_hu, double *src_hu, int64_t n0_src_hv, double *src_hv, int64_t n0_src_hB1, double *src_hB1, int64_t n0_src_hB2, double *src_hB2, int64_t n0_src_PSI, double *src_PSI, int64_t n0_src_Z, double *src_Z, int64_t n0_corio_hu, double *corio_hu, int64_t n0_corio_hv, double *corio_hv, double dtime, int64_t n0_vol, double *vol, int64_t GLM, double cpsi);
double bind_c_cpsi_global(int64_t n0_h_c, double *h_c, int64_t n0_hu_c, double *hu_c, int64_t n0_hv_c, double *hv_c, int64_t n0_hB1_c, double *hB1_c, int64_t n0_hB2_c, double *hB2_c, double cfl, int64_t n0_normal, int64_t n1_normal, double *normal, int64_t n0_mesure, double *mesure, int64_t n0_volume, double *volume, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid);
void bind_c_initialisation_swmhd(int64_t n0_h, double *h, int64_t n0_hu, double *hu, int64_t n0_hv, double *hv, int64_t n0_hB1, double *hB1, int64_t n0_hB2, double *hB2, int64_t n0_PSI, double *PSI, int64_t n0_Z, double *Z, int64_t n0_center, int64_t n1_center, double *center, int64_t choix, double k1, double k2, double eps, double tol);
void bind_c_term_source_srnh_swmhd(int64_t n0_src_h, double *src_h, int64_t n0_src_hu, double *src_hu, int64_t n0_src_hv, double *src_hv, int64_t n0_src_hB1, double *src_hB1, int64_t n0_src_hB2, double *src_hB2, int64_t n0_src_PSI, double *src_PSI, int64_t n0_src_Z, double *src_Z, int64_t n0_h_c, double *h_c, int64_t n0_hu_c, double *hu_c, int64_t n0_hv_c, double *hv_c, int64_t n0_Z_c, double *Z_c, int64_t n0_h_ghost, double *h_ghost, int64_t n0_hu_ghost, double *hu_ghost, int64_t n0_hv_ghost, double *hv_ghost, int64_t n0_Z_ghost, double *Z_ghost, int64_t n0_h_halo, double *h_halo, int64_t n0_hu_halo, double *hu_halo, int64_t n0_hv_halo, double *hv_halo, int64_t n0_Z_halo, double *Z_halo, int64_t n0_h_x, double *h_x, int64_t n0_h_y, double *h_y, int64_t n0_psi, double *psi, int64_t n0_hx_halo, double *hx_halo, int64_t n0_hy_halo, double *hy_halo, int64_t n0_psi_halo, double *psi_halo, int64_t n0_nodeidc, int64_t n1_nodeidc, int64_t *nodeidc, int64_t n0_faceidc, int64_t n1_faceidc, int64_t *faceidc, int64_t n0_cellidc, int64_t n1_cellidc, int64_t *cellidc, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_normalc, int64_t n1_normalc, int64_t n2_normalc, double *normalc, int64_t n0_namef, int64_t *namef, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_halofid, int64_t *halofid, int64_t order);
void bind_c_lf_scheme_mhd(double hu_l, double hu_r, double hv_l, double hv_r, double h_l, double h_r, double hB1_l, double hB1_r, double hB2_l, double hB2_r, double Z_l, double Z_r, int64_t n0_normal, double *normal, double mesure, double grav, int64_t n0_flux, double *flux);
void bind_c_srnh_scheme_mhd(double hu_l, double hu_r, double hv_l, double hv_r, double h_l, double h_r, double hB1_l, double hB1_r, double hB2_l, double hB2_r, double hPSI_l, double hPSI_r, double hB1c_l, double hB1c_r, double hB2c_l, double hB2c_r, double Z_l, double Z_r, int64_t n0_normal, double *normal, double mesure, double grav, int64_t n0_flux, double *flux, double cpsi);
void bind_c_explicitscheme_convective_swmhd(int64_t n0_rez_h, double *rez_h, int64_t n0_rez_hu, double *rez_hu, int64_t n0_rez_hv, double *rez_hv, int64_t n0_rez_hB1, double *rez_hB1, int64_t n0_rez_hB2, double *rez_hB2, int64_t n0_rez_PSI, double *rez_PSI, int64_t n0_rez_Z, double *rez_Z, int64_t n0_h_c, double *h_c, int64_t n0_hu_c, double *hu_c, int64_t n0_hv_c, double *hv_c, int64_t n0_hB1_c, double *hB1_c, int64_t n0_hB2_c, double *hB2_c, int64_t n0_hPSIc, double *hPSIc, int64_t n0_Z_c, double *Z_c, int64_t n0_h_ghost, double *h_ghost, int64_t n0_hu_ghost, double *hu_ghost, int64_t n0_hv_ghost, double *hv_ghost, int64_t n0_hB1_ghost, double *hB1_ghost, int64_t n0_hB2_ghost, double *hB2_ghost, int64_t n0_hPSIghost, double *hPSIghost, int64_t n0_Z_ghost, double *Z_ghost, int64_t n0_h_halo, double *h_halo, int64_t n0_hu_halo, double *hu_halo, int64_t n0_hv_halo, double *hv_halo, int64_t n0_hB1_halo, double *hB1_halo, int64_t n0_hB2_halo, double *hB2_halo, int64_t n0_hPSIhalo, double *hPSIhalo, int64_t n0_Z_halo, double *Z_halo, int64_t n0_h_x, double *h_x, int64_t n0_h_y, double *h_y, int64_t n0_hx_halo, double *hx_halo, int64_t n0_hy_halo, double *hy_halo, int64_t n0_psi, double *psi, int64_t n0_psi_halo, double *psi_halo, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_centerg, int64_t n1_centerg, double *centerg, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_mesuref, double *mesuref, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_halofid, int64_t *halofid, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_boundaryfaces, int64_t *boundaryfaces, int64_t n0_periodicboundaryfaces, int64_t *periodicboundaryfaces, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t order, double cpsi, int64_t n0_hB1_cst, double *hB1_cst, int64_t n0_hB2_cst, double *hB2_cst);

/*........................................*/


/*........................................*/

/*........................................*/
PyObject *Total_Energy_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray h_c = {.shape = NULL};
    t_ndarray hu_c = {.shape = NULL};
    t_ndarray hv_c = {.shape = NULL};
    t_ndarray hB1_c = {.shape = NULL};
    t_ndarray hB2_c = {.shape = NULL};
    t_ndarray Z_c = {.shape = NULL};
    double grav;
    t_ndarray volumec = {.shape = NULL};
    double numt;
    double numc;
    double numm;
    double nump;
    PyArrayObject *h_c_tmp;
    PyArrayObject *hu_c_tmp;
    PyArrayObject *hv_c_tmp;
    PyArrayObject *hB1_c_tmp;
    PyArrayObject *hB2_c_tmp;
    PyArrayObject *Z_c_tmp;
    PyObject *grav_tmp;
    PyArrayObject *volumec_tmp;
    PyObject *numt_tmp;
    PyObject *numc_tmp;
    PyObject *numm_tmp;
    PyObject *nump_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "h_c",
        "hu_c",
        "hv_c",
        "hB1_c",
        "hB2_c",
        "Z_c",
        "grav",
        "volumec",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!OO!", kwlist, &PyArray_Type, &h_c_tmp, &PyArray_Type, &hu_c_tmp, &PyArray_Type, &hv_c_tmp, &PyArray_Type, &hB1_c_tmp, &PyArray_Type, &hB2_c_tmp, &PyArray_Type, &Z_c_tmp, &grav_tmp, &PyArray_Type, &volumec_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(hu_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_c = pyarray_to_ndarray(hu_c_tmp);
    }
    if (!pyarray_check(hv_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_c = pyarray_to_ndarray(hv_c_tmp);
    }
    if (!pyarray_check(hB1_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_c = pyarray_to_ndarray(hB1_c_tmp);
    }
    if (!pyarray_check(hB2_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_c = pyarray_to_ndarray(hB2_c_tmp);
    }
    if (!pyarray_check(Z_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_c = pyarray_to_ndarray(Z_c_tmp);
    }
    if (PyIs_NativeFloat(grav_tmp))
    {
        grav = PyDouble_to_Double(grav_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(volumec_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volumec = pyarray_to_ndarray(volumec_tmp);
    }
    bind_c_total_energy(nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&hu_c, 0), nd_data(&hu_c), nd_ndim(&hv_c, 0), nd_data(&hv_c), nd_ndim(&hB1_c, 0), nd_data(&hB1_c), nd_ndim(&hB2_c, 0), nd_data(&hB2_c), nd_ndim(&Z_c, 0), nd_data(&Z_c), grav, nd_ndim(&volumec, 0), nd_data(&volumec), &numt, &numc, &numm, &nump);
    numt_tmp = Double_to_NumpyDouble(&numt);
    numc_tmp = Double_to_NumpyDouble(&numc);
    numm_tmp = Double_to_NumpyDouble(&numm);
    nump_tmp = Double_to_NumpyDouble(&nump);
    result = Py_BuildValue("OOOO", numt_tmp, numc_tmp, numm_tmp, nump_tmp);
    Py_DECREF(numt_tmp);
    Py_DECREF(numc_tmp);
    Py_DECREF(numm_tmp);
    Py_DECREF(nump_tmp);
    free_pointer(h_c);
    free_pointer(hu_c);
    free_pointer(hv_c);
    free_pointer(hB1_c);
    free_pointer(hB2_c);
    free_pointer(Z_c);
    free_pointer(volumec);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *explicitscheme_stabilizater_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray wx_face = {.shape = NULL};
    t_ndarray wy_face = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    double vepsilon;
    t_ndarray dissip_w = {.shape = NULL};
    t_ndarray h_c = {.shape = NULL};
    PyArrayObject *wx_face_tmp;
    PyArrayObject *wy_face_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *namef_tmp;
    PyObject *vepsilon_tmp;
    PyArrayObject *dissip_w_tmp;
    PyArrayObject *h_c_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "wx_face",
        "wy_face",
        "cellidf",
        "normalf",
        "namef",
        "vepsilon",
        "dissip_w",
        "h_c",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!OO!O!", kwlist, &PyArray_Type, &wx_face_tmp, &PyArray_Type, &wy_face_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &namef_tmp, &vepsilon_tmp, &PyArray_Type, &dissip_w_tmp, &PyArray_Type, &h_c_tmp))
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
    if (PyIs_NativeFloat(vepsilon_tmp))
    {
        vepsilon = PyDouble_to_Double(vepsilon_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(dissip_w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dissip_w = pyarray_to_ndarray(dissip_w_tmp);
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    bind_c_explicitscheme_stabilizater(nd_ndim(&wx_face, 0), nd_data(&wx_face), nd_ndim(&wy_face, 0), nd_data(&wy_face), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&namef, 0), nd_data(&namef), vepsilon, nd_ndim(&dissip_w, 0), nd_data(&dissip_w), nd_ndim(&h_c, 0), nd_data(&h_c));
    result = Py_BuildValue("");
    free_pointer(wx_face);
    free_pointer(wy_face);
    free_pointer(cellidf);
    free_pointer(normalf);
    free_pointer(namef);
    free_pointer(dissip_w);
    free_pointer(h_c);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *time_step_SWMHD_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray h_c = {.shape = NULL};
    t_ndarray hu_c = {.shape = NULL};
    t_ndarray hv_c = {.shape = NULL};
    t_ndarray hB1_c = {.shape = NULL};
    t_ndarray hB2_c = {.shape = NULL};
    double cfl;
    t_ndarray normal = {.shape = NULL};
    t_ndarray mesure = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray faceid = {.shape = NULL};
    double dt;
    PyArrayObject *h_c_tmp;
    PyArrayObject *hu_c_tmp;
    PyArrayObject *hv_c_tmp;
    PyArrayObject *hB1_c_tmp;
    PyArrayObject *hB2_c_tmp;
    PyObject *cfl_tmp;
    PyArrayObject *normal_tmp;
    PyArrayObject *mesure_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *faceid_tmp;
    PyObject *dt_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "h_c",
        "hu_c",
        "hv_c",
        "hB1_c",
        "hB2_c",
        "cfl",
        "normal",
        "mesure",
        "volume",
        "faceid",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!OO!O!O!O!", kwlist, &PyArray_Type, &h_c_tmp, &PyArray_Type, &hu_c_tmp, &PyArray_Type, &hv_c_tmp, &PyArray_Type, &hB1_c_tmp, &PyArray_Type, &hB2_c_tmp, &cfl_tmp, &PyArray_Type, &normal_tmp, &PyArray_Type, &mesure_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &faceid_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(hu_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_c = pyarray_to_ndarray(hu_c_tmp);
    }
    if (!pyarray_check(hv_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_c = pyarray_to_ndarray(hv_c_tmp);
    }
    if (!pyarray_check(hB1_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_c = pyarray_to_ndarray(hB1_c_tmp);
    }
    if (!pyarray_check(hB2_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_c = pyarray_to_ndarray(hB2_c_tmp);
    }
    if (PyIs_NativeFloat(cfl_tmp))
    {
        cfl = PyDouble_to_Double(cfl_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(normal_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normal = pyarray_to_ndarray(normal_tmp);
    }
    if (!pyarray_check(mesure_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesure = pyarray_to_ndarray(mesure_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = pyarray_to_ndarray(faceid_tmp);
    }
    dt = bind_c_time_step_swmhd(nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&hu_c, 0), nd_data(&hu_c), nd_ndim(&hv_c, 0), nd_data(&hv_c), nd_ndim(&hB1_c, 0), nd_data(&hB1_c), nd_ndim(&hB2_c, 0), nd_data(&hB2_c), cfl, nd_ndim(&normal, 0), nd_ndim(&normal, 1), nd_data(&normal), nd_ndim(&mesure, 0), nd_data(&mesure), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid));
    dt_tmp = Double_to_NumpyDouble(&dt);
    result = Py_BuildValue("O", dt_tmp);
    Py_DECREF(dt_tmp);
    free_pointer(h_c);
    free_pointer(hu_c);
    free_pointer(hv_c);
    free_pointer(hB1_c);
    free_pointer(hB2_c);
    free_pointer(normal);
    free_pointer(mesure);
    free_pointer(volume);
    free_pointer(faceid);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *update_SWMHD_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray h_c = {.shape = NULL};
    t_ndarray hu_c = {.shape = NULL};
    t_ndarray hv_c = {.shape = NULL};
    t_ndarray hB1_c = {.shape = NULL};
    t_ndarray hB2_c = {.shape = NULL};
    t_ndarray PSI_c = {.shape = NULL};
    t_ndarray Z_c = {.shape = NULL};
    t_ndarray rez_h = {.shape = NULL};
    t_ndarray rez_hu = {.shape = NULL};
    t_ndarray rez_hv = {.shape = NULL};
    t_ndarray rez_hB1 = {.shape = NULL};
    t_ndarray rez_hB2 = {.shape = NULL};
    t_ndarray rez_PSI = {.shape = NULL};
    t_ndarray rez_Z = {.shape = NULL};
    t_ndarray src_h = {.shape = NULL};
    t_ndarray src_hu = {.shape = NULL};
    t_ndarray src_hv = {.shape = NULL};
    t_ndarray src_hB1 = {.shape = NULL};
    t_ndarray src_hB2 = {.shape = NULL};
    t_ndarray src_PSI = {.shape = NULL};
    t_ndarray src_Z = {.shape = NULL};
    t_ndarray corio_hu = {.shape = NULL};
    t_ndarray corio_hv = {.shape = NULL};
    double dtime;
    t_ndarray vol = {.shape = NULL};
    int64_t GLM;
    double cpsi;
    PyArrayObject *h_c_tmp;
    PyArrayObject *hu_c_tmp;
    PyArrayObject *hv_c_tmp;
    PyArrayObject *hB1_c_tmp;
    PyArrayObject *hB2_c_tmp;
    PyArrayObject *PSI_c_tmp;
    PyArrayObject *Z_c_tmp;
    PyArrayObject *rez_h_tmp;
    PyArrayObject *rez_hu_tmp;
    PyArrayObject *rez_hv_tmp;
    PyArrayObject *rez_hB1_tmp;
    PyArrayObject *rez_hB2_tmp;
    PyArrayObject *rez_PSI_tmp;
    PyArrayObject *rez_Z_tmp;
    PyArrayObject *src_h_tmp;
    PyArrayObject *src_hu_tmp;
    PyArrayObject *src_hv_tmp;
    PyArrayObject *src_hB1_tmp;
    PyArrayObject *src_hB2_tmp;
    PyArrayObject *src_PSI_tmp;
    PyArrayObject *src_Z_tmp;
    PyArrayObject *corio_hu_tmp;
    PyArrayObject *corio_hv_tmp;
    PyObject *dtime_tmp;
    PyArrayObject *vol_tmp;
    PyObject *GLM_tmp;
    PyObject *cpsi_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "h_c",
        "hu_c",
        "hv_c",
        "hB1_c",
        "hB2_c",
        "PSI_c",
        "Z_c",
        "rez_h",
        "rez_hu",
        "rez_hv",
        "rez_hB1",
        "rez_hB2",
        "rez_PSI",
        "rez_Z",
        "src_h",
        "src_hu",
        "src_hv",
        "src_hB1",
        "src_hB2",
        "src_PSI",
        "src_Z",
        "corio_hu",
        "corio_hv",
        "dtime",
        "vol",
        "GLM",
        "cpsi",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!OO!OO", kwlist, &PyArray_Type, &h_c_tmp, &PyArray_Type, &hu_c_tmp, &PyArray_Type, &hv_c_tmp, &PyArray_Type, &hB1_c_tmp, &PyArray_Type, &hB2_c_tmp, &PyArray_Type, &PSI_c_tmp, &PyArray_Type, &Z_c_tmp, &PyArray_Type, &rez_h_tmp, &PyArray_Type, &rez_hu_tmp, &PyArray_Type, &rez_hv_tmp, &PyArray_Type, &rez_hB1_tmp, &PyArray_Type, &rez_hB2_tmp, &PyArray_Type, &rez_PSI_tmp, &PyArray_Type, &rez_Z_tmp, &PyArray_Type, &src_h_tmp, &PyArray_Type, &src_hu_tmp, &PyArray_Type, &src_hv_tmp, &PyArray_Type, &src_hB1_tmp, &PyArray_Type, &src_hB2_tmp, &PyArray_Type, &src_PSI_tmp, &PyArray_Type, &src_Z_tmp, &PyArray_Type, &corio_hu_tmp, &PyArray_Type, &corio_hv_tmp, &dtime_tmp, &PyArray_Type, &vol_tmp, &GLM_tmp, &cpsi_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(hu_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_c = pyarray_to_ndarray(hu_c_tmp);
    }
    if (!pyarray_check(hv_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_c = pyarray_to_ndarray(hv_c_tmp);
    }
    if (!pyarray_check(hB1_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_c = pyarray_to_ndarray(hB1_c_tmp);
    }
    if (!pyarray_check(hB2_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_c = pyarray_to_ndarray(hB2_c_tmp);
    }
    if (!pyarray_check(PSI_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        PSI_c = pyarray_to_ndarray(PSI_c_tmp);
    }
    if (!pyarray_check(Z_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_c = pyarray_to_ndarray(Z_c_tmp);
    }
    if (!pyarray_check(rez_h_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_h = pyarray_to_ndarray(rez_h_tmp);
    }
    if (!pyarray_check(rez_hu_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hu = pyarray_to_ndarray(rez_hu_tmp);
    }
    if (!pyarray_check(rez_hv_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hv = pyarray_to_ndarray(rez_hv_tmp);
    }
    if (!pyarray_check(rez_hB1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hB1 = pyarray_to_ndarray(rez_hB1_tmp);
    }
    if (!pyarray_check(rez_hB2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hB2 = pyarray_to_ndarray(rez_hB2_tmp);
    }
    if (!pyarray_check(rez_PSI_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_PSI = pyarray_to_ndarray(rez_PSI_tmp);
    }
    if (!pyarray_check(rez_Z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_Z = pyarray_to_ndarray(rez_Z_tmp);
    }
    if (!pyarray_check(src_h_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_h = pyarray_to_ndarray(src_h_tmp);
    }
    if (!pyarray_check(src_hu_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hu = pyarray_to_ndarray(src_hu_tmp);
    }
    if (!pyarray_check(src_hv_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hv = pyarray_to_ndarray(src_hv_tmp);
    }
    if (!pyarray_check(src_hB1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hB1 = pyarray_to_ndarray(src_hB1_tmp);
    }
    if (!pyarray_check(src_hB2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hB2 = pyarray_to_ndarray(src_hB2_tmp);
    }
    if (!pyarray_check(src_PSI_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_PSI = pyarray_to_ndarray(src_PSI_tmp);
    }
    if (!pyarray_check(src_Z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_Z = pyarray_to_ndarray(src_Z_tmp);
    }
    if (!pyarray_check(corio_hu_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        corio_hu = pyarray_to_ndarray(corio_hu_tmp);
    }
    if (!pyarray_check(corio_hv_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        corio_hv = pyarray_to_ndarray(corio_hv_tmp);
    }
    if (PyIs_NativeFloat(dtime_tmp))
    {
        dtime = PyDouble_to_Double(dtime_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(vol_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        vol = pyarray_to_ndarray(vol_tmp);
    }
    if (PyIs_NativeInt(GLM_tmp))
    {
        GLM = PyInt64_to_Int64(GLM_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (PyIs_NativeFloat(cpsi_tmp))
    {
        cpsi = PyDouble_to_Double(cpsi_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    bind_c_update_swmhd(nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&hu_c, 0), nd_data(&hu_c), nd_ndim(&hv_c, 0), nd_data(&hv_c), nd_ndim(&hB1_c, 0), nd_data(&hB1_c), nd_ndim(&hB2_c, 0), nd_data(&hB2_c), nd_ndim(&PSI_c, 0), nd_data(&PSI_c), nd_ndim(&Z_c, 0), nd_data(&Z_c), nd_ndim(&rez_h, 0), nd_data(&rez_h), nd_ndim(&rez_hu, 0), nd_data(&rez_hu), nd_ndim(&rez_hv, 0), nd_data(&rez_hv), nd_ndim(&rez_hB1, 0), nd_data(&rez_hB1), nd_ndim(&rez_hB2, 0), nd_data(&rez_hB2), nd_ndim(&rez_PSI, 0), nd_data(&rez_PSI), nd_ndim(&rez_Z, 0), nd_data(&rez_Z), nd_ndim(&src_h, 0), nd_data(&src_h), nd_ndim(&src_hu, 0), nd_data(&src_hu), nd_ndim(&src_hv, 0), nd_data(&src_hv), nd_ndim(&src_hB1, 0), nd_data(&src_hB1), nd_ndim(&src_hB2, 0), nd_data(&src_hB2), nd_ndim(&src_PSI, 0), nd_data(&src_PSI), nd_ndim(&src_Z, 0), nd_data(&src_Z), nd_ndim(&corio_hu, 0), nd_data(&corio_hu), nd_ndim(&corio_hv, 0), nd_data(&corio_hv), dtime, nd_ndim(&vol, 0), nd_data(&vol), GLM, cpsi);
    result = Py_BuildValue("");
    free_pointer(h_c);
    free_pointer(hu_c);
    free_pointer(hv_c);
    free_pointer(hB1_c);
    free_pointer(hB2_c);
    free_pointer(PSI_c);
    free_pointer(Z_c);
    free_pointer(rez_h);
    free_pointer(rez_hu);
    free_pointer(rez_hv);
    free_pointer(rez_hB1);
    free_pointer(rez_hB2);
    free_pointer(rez_PSI);
    free_pointer(rez_Z);
    free_pointer(src_h);
    free_pointer(src_hu);
    free_pointer(src_hv);
    free_pointer(src_hB1);
    free_pointer(src_hB2);
    free_pointer(src_PSI);
    free_pointer(src_Z);
    free_pointer(corio_hu);
    free_pointer(corio_hv);
    free_pointer(vol);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *cpsi_global_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray h_c = {.shape = NULL};
    t_ndarray hu_c = {.shape = NULL};
    t_ndarray hv_c = {.shape = NULL};
    t_ndarray hB1_c = {.shape = NULL};
    t_ndarray hB2_c = {.shape = NULL};
    double cfl;
    t_ndarray normal = {.shape = NULL};
    t_ndarray mesure = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray faceid = {.shape = NULL};
    double cpsiglobal;
    PyArrayObject *h_c_tmp;
    PyArrayObject *hu_c_tmp;
    PyArrayObject *hv_c_tmp;
    PyArrayObject *hB1_c_tmp;
    PyArrayObject *hB2_c_tmp;
    PyObject *cfl_tmp;
    PyArrayObject *normal_tmp;
    PyArrayObject *mesure_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *faceid_tmp;
    PyObject *cpsiglobal_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "h_c",
        "hu_c",
        "hv_c",
        "hB1_c",
        "hB2_c",
        "cfl",
        "normal",
        "mesure",
        "volume",
        "faceid",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!OO!O!O!O!", kwlist, &PyArray_Type, &h_c_tmp, &PyArray_Type, &hu_c_tmp, &PyArray_Type, &hv_c_tmp, &PyArray_Type, &hB1_c_tmp, &PyArray_Type, &hB2_c_tmp, &cfl_tmp, &PyArray_Type, &normal_tmp, &PyArray_Type, &mesure_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &faceid_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(hu_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_c = pyarray_to_ndarray(hu_c_tmp);
    }
    if (!pyarray_check(hv_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_c = pyarray_to_ndarray(hv_c_tmp);
    }
    if (!pyarray_check(hB1_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_c = pyarray_to_ndarray(hB1_c_tmp);
    }
    if (!pyarray_check(hB2_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_c = pyarray_to_ndarray(hB2_c_tmp);
    }
    if (PyIs_NativeFloat(cfl_tmp))
    {
        cfl = PyDouble_to_Double(cfl_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(normal_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normal = pyarray_to_ndarray(normal_tmp);
    }
    if (!pyarray_check(mesure_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesure = pyarray_to_ndarray(mesure_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = pyarray_to_ndarray(faceid_tmp);
    }
    cpsiglobal = bind_c_cpsi_global(nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&hu_c, 0), nd_data(&hu_c), nd_ndim(&hv_c, 0), nd_data(&hv_c), nd_ndim(&hB1_c, 0), nd_data(&hB1_c), nd_ndim(&hB2_c, 0), nd_data(&hB2_c), cfl, nd_ndim(&normal, 0), nd_ndim(&normal, 1), nd_data(&normal), nd_ndim(&mesure, 0), nd_data(&mesure), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid));
    cpsiglobal_tmp = Double_to_NumpyDouble(&cpsiglobal);
    result = Py_BuildValue("O", cpsiglobal_tmp);
    Py_DECREF(cpsiglobal_tmp);
    free_pointer(h_c);
    free_pointer(hu_c);
    free_pointer(hv_c);
    free_pointer(hB1_c);
    free_pointer(hB2_c);
    free_pointer(normal);
    free_pointer(mesure);
    free_pointer(volume);
    free_pointer(faceid);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *initialisation_SWMHD_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray h = {.shape = NULL};
    t_ndarray hu = {.shape = NULL};
    t_ndarray hv = {.shape = NULL};
    t_ndarray hB1 = {.shape = NULL};
    t_ndarray hB2 = {.shape = NULL};
    t_ndarray PSI = {.shape = NULL};
    t_ndarray Z = {.shape = NULL};
    t_ndarray center = {.shape = NULL};
    int64_t choix;
    double k1;
    double k2;
    double eps;
    double tol;
    PyArrayObject *h_tmp;
    PyArrayObject *hu_tmp;
    PyArrayObject *hv_tmp;
    PyArrayObject *hB1_tmp;
    PyArrayObject *hB2_tmp;
    PyArrayObject *PSI_tmp;
    PyArrayObject *Z_tmp;
    PyArrayObject *center_tmp;
    PyObject *choix_tmp;
    PyObject *k1_tmp;
    PyObject *k2_tmp;
    PyObject *eps_tmp;
    PyObject *tol_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "h",
        "hu",
        "hv",
        "hB1",
        "hB2",
        "PSI",
        "Z",
        "center",
        "choix",
        "k1",
        "k2",
        "eps",
        "tol",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!OOOOO", kwlist, &PyArray_Type, &h_tmp, &PyArray_Type, &hu_tmp, &PyArray_Type, &hv_tmp, &PyArray_Type, &hB1_tmp, &PyArray_Type, &hB2_tmp, &PyArray_Type, &PSI_tmp, &PyArray_Type, &Z_tmp, &PyArray_Type, &center_tmp, &choix_tmp, &k1_tmp, &k2_tmp, &eps_tmp, &tol_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(h_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h = pyarray_to_ndarray(h_tmp);
    }
    if (!pyarray_check(hu_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu = pyarray_to_ndarray(hu_tmp);
    }
    if (!pyarray_check(hv_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv = pyarray_to_ndarray(hv_tmp);
    }
    if (!pyarray_check(hB1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1 = pyarray_to_ndarray(hB1_tmp);
    }
    if (!pyarray_check(hB2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2 = pyarray_to_ndarray(hB2_tmp);
    }
    if (!pyarray_check(PSI_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        PSI = pyarray_to_ndarray(PSI_tmp);
    }
    if (!pyarray_check(Z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z = pyarray_to_ndarray(Z_tmp);
    }
    if (!pyarray_check(center_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        center = pyarray_to_ndarray(center_tmp);
    }
    if (PyIs_NativeInt(choix_tmp))
    {
        choix = PyInt64_to_Int64(choix_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (PyIs_NativeFloat(k1_tmp))
    {
        k1 = PyDouble_to_Double(k1_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(k2_tmp))
    {
        k2 = PyDouble_to_Double(k2_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(eps_tmp))
    {
        eps = PyDouble_to_Double(eps_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(tol_tmp))
    {
        tol = PyDouble_to_Double(tol_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    bind_c_initialisation_swmhd(nd_ndim(&h, 0), nd_data(&h), nd_ndim(&hu, 0), nd_data(&hu), nd_ndim(&hv, 0), nd_data(&hv), nd_ndim(&hB1, 0), nd_data(&hB1), nd_ndim(&hB2, 0), nd_data(&hB2), nd_ndim(&PSI, 0), nd_data(&PSI), nd_ndim(&Z, 0), nd_data(&Z), nd_ndim(&center, 0), nd_ndim(&center, 1), nd_data(&center), choix, k1, k2, eps, tol);
    result = Py_BuildValue("");
    free_pointer(h);
    free_pointer(hu);
    free_pointer(hv);
    free_pointer(hB1);
    free_pointer(hB2);
    free_pointer(PSI);
    free_pointer(Z);
    free_pointer(center);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *term_source_srnh_SWMHD_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray src_h = {.shape = NULL};
    t_ndarray src_hu = {.shape = NULL};
    t_ndarray src_hv = {.shape = NULL};
    t_ndarray src_hB1 = {.shape = NULL};
    t_ndarray src_hB2 = {.shape = NULL};
    t_ndarray src_PSI = {.shape = NULL};
    t_ndarray src_Z = {.shape = NULL};
    t_ndarray h_c = {.shape = NULL};
    t_ndarray hu_c = {.shape = NULL};
    t_ndarray hv_c = {.shape = NULL};
    t_ndarray Z_c = {.shape = NULL};
    t_ndarray h_ghost = {.shape = NULL};
    t_ndarray hu_ghost = {.shape = NULL};
    t_ndarray hv_ghost = {.shape = NULL};
    t_ndarray Z_ghost = {.shape = NULL};
    t_ndarray h_halo = {.shape = NULL};
    t_ndarray hu_halo = {.shape = NULL};
    t_ndarray hv_halo = {.shape = NULL};
    t_ndarray Z_halo = {.shape = NULL};
    t_ndarray h_x = {.shape = NULL};
    t_ndarray h_y = {.shape = NULL};
    t_ndarray psi = {.shape = NULL};
    t_ndarray hx_halo = {.shape = NULL};
    t_ndarray hy_halo = {.shape = NULL};
    t_ndarray psi_halo = {.shape = NULL};
    t_ndarray nodeidc = {.shape = NULL};
    t_ndarray faceidc = {.shape = NULL};
    t_ndarray cellidc = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray normalc = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    int64_t order;
    PyArrayObject *src_h_tmp;
    PyArrayObject *src_hu_tmp;
    PyArrayObject *src_hv_tmp;
    PyArrayObject *src_hB1_tmp;
    PyArrayObject *src_hB2_tmp;
    PyArrayObject *src_PSI_tmp;
    PyArrayObject *src_Z_tmp;
    PyArrayObject *h_c_tmp;
    PyArrayObject *hu_c_tmp;
    PyArrayObject *hv_c_tmp;
    PyArrayObject *Z_c_tmp;
    PyArrayObject *h_ghost_tmp;
    PyArrayObject *hu_ghost_tmp;
    PyArrayObject *hv_ghost_tmp;
    PyArrayObject *Z_ghost_tmp;
    PyArrayObject *h_halo_tmp;
    PyArrayObject *hu_halo_tmp;
    PyArrayObject *hv_halo_tmp;
    PyArrayObject *Z_halo_tmp;
    PyArrayObject *h_x_tmp;
    PyArrayObject *h_y_tmp;
    PyArrayObject *psi_tmp;
    PyArrayObject *hx_halo_tmp;
    PyArrayObject *hy_halo_tmp;
    PyArrayObject *psi_halo_tmp;
    PyArrayObject *nodeidc_tmp;
    PyArrayObject *faceidc_tmp;
    PyArrayObject *cellidc_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *normalc_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *centerf_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *halofid_tmp;
    PyObject *order_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "src_h",
        "src_hu",
        "src_hv",
        "src_hB1",
        "src_hB2",
        "src_PSI",
        "src_Z",
        "h_c",
        "hu_c",
        "hv_c",
        "Z_c",
        "h_ghost",
        "hu_ghost",
        "hv_ghost",
        "Z_ghost",
        "h_halo",
        "hu_halo",
        "hv_halo",
        "Z_halo",
        "h_x",
        "h_y",
        "psi",
        "hx_halo",
        "hy_halo",
        "psi_halo",
        "nodeidc",
        "faceidc",
        "cellidc",
        "cellidf",
        "centerc",
        "normalc",
        "namef",
        "centerf",
        "centerh",
        "vertexn",
        "halofid",
        "order",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O", kwlist, &PyArray_Type, &src_h_tmp, &PyArray_Type, &src_hu_tmp, &PyArray_Type, &src_hv_tmp, &PyArray_Type, &src_hB1_tmp, &PyArray_Type, &src_hB2_tmp, &PyArray_Type, &src_PSI_tmp, &PyArray_Type, &src_Z_tmp, &PyArray_Type, &h_c_tmp, &PyArray_Type, &hu_c_tmp, &PyArray_Type, &hv_c_tmp, &PyArray_Type, &Z_c_tmp, &PyArray_Type, &h_ghost_tmp, &PyArray_Type, &hu_ghost_tmp, &PyArray_Type, &hv_ghost_tmp, &PyArray_Type, &Z_ghost_tmp, &PyArray_Type, &h_halo_tmp, &PyArray_Type, &hu_halo_tmp, &PyArray_Type, &hv_halo_tmp, &PyArray_Type, &Z_halo_tmp, &PyArray_Type, &h_x_tmp, &PyArray_Type, &h_y_tmp, &PyArray_Type, &psi_tmp, &PyArray_Type, &hx_halo_tmp, &PyArray_Type, &hy_halo_tmp, &PyArray_Type, &psi_halo_tmp, &PyArray_Type, &nodeidc_tmp, &PyArray_Type, &faceidc_tmp, &PyArray_Type, &cellidc_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &normalc_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &centerf_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &halofid_tmp, &order_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(src_h_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_h = pyarray_to_ndarray(src_h_tmp);
    }
    if (!pyarray_check(src_hu_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hu = pyarray_to_ndarray(src_hu_tmp);
    }
    if (!pyarray_check(src_hv_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hv = pyarray_to_ndarray(src_hv_tmp);
    }
    if (!pyarray_check(src_hB1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hB1 = pyarray_to_ndarray(src_hB1_tmp);
    }
    if (!pyarray_check(src_hB2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_hB2 = pyarray_to_ndarray(src_hB2_tmp);
    }
    if (!pyarray_check(src_PSI_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_PSI = pyarray_to_ndarray(src_PSI_tmp);
    }
    if (!pyarray_check(src_Z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_Z = pyarray_to_ndarray(src_Z_tmp);
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(hu_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_c = pyarray_to_ndarray(hu_c_tmp);
    }
    if (!pyarray_check(hv_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_c = pyarray_to_ndarray(hv_c_tmp);
    }
    if (!pyarray_check(Z_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_c = pyarray_to_ndarray(Z_c_tmp);
    }
    if (!pyarray_check(h_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_ghost = pyarray_to_ndarray(h_ghost_tmp);
    }
    if (!pyarray_check(hu_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_ghost = pyarray_to_ndarray(hu_ghost_tmp);
    }
    if (!pyarray_check(hv_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_ghost = pyarray_to_ndarray(hv_ghost_tmp);
    }
    if (!pyarray_check(Z_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_ghost = pyarray_to_ndarray(Z_ghost_tmp);
    }
    if (!pyarray_check(h_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_halo = pyarray_to_ndarray(h_halo_tmp);
    }
    if (!pyarray_check(hu_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_halo = pyarray_to_ndarray(hu_halo_tmp);
    }
    if (!pyarray_check(hv_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_halo = pyarray_to_ndarray(hv_halo_tmp);
    }
    if (!pyarray_check(Z_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_halo = pyarray_to_ndarray(Z_halo_tmp);
    }
    if (!pyarray_check(h_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_x = pyarray_to_ndarray(h_x_tmp);
    }
    if (!pyarray_check(h_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_y = pyarray_to_ndarray(h_y_tmp);
    }
    if (!pyarray_check(psi_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi = pyarray_to_ndarray(psi_tmp);
    }
    if (!pyarray_check(hx_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hx_halo = pyarray_to_ndarray(hx_halo_tmp);
    }
    if (!pyarray_check(hy_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hy_halo = pyarray_to_ndarray(hy_halo_tmp);
    }
    if (!pyarray_check(psi_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi_halo = pyarray_to_ndarray(psi_halo_tmp);
    }
    if (!pyarray_check(nodeidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidc = pyarray_to_ndarray(nodeidc_tmp);
    }
    if (!pyarray_check(faceidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceidc = pyarray_to_ndarray(faceidc_tmp);
    }
    if (!pyarray_check(cellidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidc = pyarray_to_ndarray(cellidc_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(normalc_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalc = pyarray_to_ndarray(normalc_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
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
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
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
    bind_c_term_source_srnh_swmhd(nd_ndim(&src_h, 0), nd_data(&src_h), nd_ndim(&src_hu, 0), nd_data(&src_hu), nd_ndim(&src_hv, 0), nd_data(&src_hv), nd_ndim(&src_hB1, 0), nd_data(&src_hB1), nd_ndim(&src_hB2, 0), nd_data(&src_hB2), nd_ndim(&src_PSI, 0), nd_data(&src_PSI), nd_ndim(&src_Z, 0), nd_data(&src_Z), nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&hu_c, 0), nd_data(&hu_c), nd_ndim(&hv_c, 0), nd_data(&hv_c), nd_ndim(&Z_c, 0), nd_data(&Z_c), nd_ndim(&h_ghost, 0), nd_data(&h_ghost), nd_ndim(&hu_ghost, 0), nd_data(&hu_ghost), nd_ndim(&hv_ghost, 0), nd_data(&hv_ghost), nd_ndim(&Z_ghost, 0), nd_data(&Z_ghost), nd_ndim(&h_halo, 0), nd_data(&h_halo), nd_ndim(&hu_halo, 0), nd_data(&hu_halo), nd_ndim(&hv_halo, 0), nd_data(&hv_halo), nd_ndim(&Z_halo, 0), nd_data(&Z_halo), nd_ndim(&h_x, 0), nd_data(&h_x), nd_ndim(&h_y, 0), nd_data(&h_y), nd_ndim(&psi, 0), nd_data(&psi), nd_ndim(&hx_halo, 0), nd_data(&hx_halo), nd_ndim(&hy_halo, 0), nd_data(&hy_halo), nd_ndim(&psi_halo, 0), nd_data(&psi_halo), nd_ndim(&nodeidc, 0), nd_ndim(&nodeidc, 1), nd_data(&nodeidc), nd_ndim(&faceidc, 0), nd_ndim(&faceidc, 1), nd_data(&faceidc), nd_ndim(&cellidc, 0), nd_ndim(&cellidc, 1), nd_data(&cellidc), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&normalc, 0), nd_ndim(&normalc, 1), nd_ndim(&normalc, 2), nd_data(&normalc), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&halofid, 0), nd_data(&halofid), order);
    result = Py_BuildValue("");
    free_pointer(src_h);
    free_pointer(src_hu);
    free_pointer(src_hv);
    free_pointer(src_hB1);
    free_pointer(src_hB2);
    free_pointer(src_PSI);
    free_pointer(src_Z);
    free_pointer(h_c);
    free_pointer(hu_c);
    free_pointer(hv_c);
    free_pointer(Z_c);
    free_pointer(h_ghost);
    free_pointer(hu_ghost);
    free_pointer(hv_ghost);
    free_pointer(Z_ghost);
    free_pointer(h_halo);
    free_pointer(hu_halo);
    free_pointer(hv_halo);
    free_pointer(Z_halo);
    free_pointer(h_x);
    free_pointer(h_y);
    free_pointer(psi);
    free_pointer(hx_halo);
    free_pointer(hy_halo);
    free_pointer(psi_halo);
    free_pointer(nodeidc);
    free_pointer(faceidc);
    free_pointer(cellidc);
    free_pointer(cellidf);
    free_pointer(centerc);
    free_pointer(normalc);
    free_pointer(namef);
    free_pointer(centerf);
    free_pointer(centerh);
    free_pointer(vertexn);
    free_pointer(halofid);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *LF_scheme_MHD_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    double hu_l;
    double hu_r;
    double hv_l;
    double hv_r;
    double h_l;
    double h_r;
    double hB1_l;
    double hB1_r;
    double hB2_l;
    double hB2_r;
    double Z_l;
    double Z_r;
    t_ndarray normal = {.shape = NULL};
    double mesure;
    double grav;
    t_ndarray flux = {.shape = NULL};
    PyObject *hu_l_tmp;
    PyObject *hu_r_tmp;
    PyObject *hv_l_tmp;
    PyObject *hv_r_tmp;
    PyObject *h_l_tmp;
    PyObject *h_r_tmp;
    PyObject *hB1_l_tmp;
    PyObject *hB1_r_tmp;
    PyObject *hB2_l_tmp;
    PyObject *hB2_r_tmp;
    PyObject *Z_l_tmp;
    PyObject *Z_r_tmp;
    PyArrayObject *normal_tmp;
    PyObject *mesure_tmp;
    PyObject *grav_tmp;
    PyArrayObject *flux_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "hu_l",
        "hu_r",
        "hv_l",
        "hv_r",
        "h_l",
        "h_r",
        "hB1_l",
        "hB1_r",
        "hB2_l",
        "hB2_r",
        "Z_l",
        "Z_r",
        "normal",
        "mesure",
        "grav",
        "flux",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOOOOOO!OOO!", kwlist, &hu_l_tmp, &hu_r_tmp, &hv_l_tmp, &hv_r_tmp, &h_l_tmp, &h_r_tmp, &hB1_l_tmp, &hB1_r_tmp, &hB2_l_tmp, &hB2_r_tmp, &Z_l_tmp, &Z_r_tmp, &PyArray_Type, &normal_tmp, &mesure_tmp, &grav_tmp, &PyArray_Type, &flux_tmp))
    {
        return NULL;
    }
    if (PyIs_NativeFloat(hu_l_tmp))
    {
        hu_l = PyDouble_to_Double(hu_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hu_r_tmp))
    {
        hu_r = PyDouble_to_Double(hu_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hv_l_tmp))
    {
        hv_l = PyDouble_to_Double(hv_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hv_r_tmp))
    {
        hv_r = PyDouble_to_Double(hv_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(h_l_tmp))
    {
        h_l = PyDouble_to_Double(h_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(h_r_tmp))
    {
        h_r = PyDouble_to_Double(h_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB1_l_tmp))
    {
        hB1_l = PyDouble_to_Double(hB1_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB1_r_tmp))
    {
        hB1_r = PyDouble_to_Double(hB1_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB2_l_tmp))
    {
        hB2_l = PyDouble_to_Double(hB2_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB2_r_tmp))
    {
        hB2_r = PyDouble_to_Double(hB2_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(Z_l_tmp))
    {
        Z_l = PyDouble_to_Double(Z_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(Z_r_tmp))
    {
        Z_r = PyDouble_to_Double(Z_r_tmp);
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
    if (PyIs_NativeFloat(mesure_tmp))
    {
        mesure = PyDouble_to_Double(mesure_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(grav_tmp))
    {
        grav = PyDouble_to_Double(grav_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(flux_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        flux = pyarray_to_ndarray(flux_tmp);
    }
    bind_c_lf_scheme_mhd(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, nd_ndim(&normal, 0), nd_data(&normal), mesure, grav, nd_ndim(&flux, 0), nd_data(&flux));
    result = Py_BuildValue("");
    free_pointer(normal);
    free_pointer(flux);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *srnh_scheme_MHD_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    double hu_l;
    double hu_r;
    double hv_l;
    double hv_r;
    double h_l;
    double h_r;
    double hB1_l;
    double hB1_r;
    double hB2_l;
    double hB2_r;
    double hPSI_l;
    double hPSI_r;
    double hB1c_l;
    double hB1c_r;
    double hB2c_l;
    double hB2c_r;
    double Z_l;
    double Z_r;
    t_ndarray normal = {.shape = NULL};
    double mesure;
    double grav;
    t_ndarray flux = {.shape = NULL};
    double cpsi;
    PyObject *hu_l_tmp;
    PyObject *hu_r_tmp;
    PyObject *hv_l_tmp;
    PyObject *hv_r_tmp;
    PyObject *h_l_tmp;
    PyObject *h_r_tmp;
    PyObject *hB1_l_tmp;
    PyObject *hB1_r_tmp;
    PyObject *hB2_l_tmp;
    PyObject *hB2_r_tmp;
    PyObject *hPSI_l_tmp;
    PyObject *hPSI_r_tmp;
    PyObject *hB1c_l_tmp;
    PyObject *hB1c_r_tmp;
    PyObject *hB2c_l_tmp;
    PyObject *hB2c_r_tmp;
    PyObject *Z_l_tmp;
    PyObject *Z_r_tmp;
    PyArrayObject *normal_tmp;
    PyObject *mesure_tmp;
    PyObject *grav_tmp;
    PyArrayObject *flux_tmp;
    PyObject *cpsi_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "hu_l",
        "hu_r",
        "hv_l",
        "hv_r",
        "h_l",
        "h_r",
        "hB1_l",
        "hB1_r",
        "hB2_l",
        "hB2_r",
        "hPSI_l",
        "hPSI_r",
        "hB1c_l",
        "hB1c_r",
        "hB2c_l",
        "hB2c_r",
        "Z_l",
        "Z_r",
        "normal",
        "mesure",
        "grav",
        "flux",
        "cpsi",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOOOOOOOOOOOO!OOO!O", kwlist, &hu_l_tmp, &hu_r_tmp, &hv_l_tmp, &hv_r_tmp, &h_l_tmp, &h_r_tmp, &hB1_l_tmp, &hB1_r_tmp, &hB2_l_tmp, &hB2_r_tmp, &hPSI_l_tmp, &hPSI_r_tmp, &hB1c_l_tmp, &hB1c_r_tmp, &hB2c_l_tmp, &hB2c_r_tmp, &Z_l_tmp, &Z_r_tmp, &PyArray_Type, &normal_tmp, &mesure_tmp, &grav_tmp, &PyArray_Type, &flux_tmp, &cpsi_tmp))
    {
        return NULL;
    }
    if (PyIs_NativeFloat(hu_l_tmp))
    {
        hu_l = PyDouble_to_Double(hu_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hu_r_tmp))
    {
        hu_r = PyDouble_to_Double(hu_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hv_l_tmp))
    {
        hv_l = PyDouble_to_Double(hv_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hv_r_tmp))
    {
        hv_r = PyDouble_to_Double(hv_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(h_l_tmp))
    {
        h_l = PyDouble_to_Double(h_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(h_r_tmp))
    {
        h_r = PyDouble_to_Double(h_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB1_l_tmp))
    {
        hB1_l = PyDouble_to_Double(hB1_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB1_r_tmp))
    {
        hB1_r = PyDouble_to_Double(hB1_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB2_l_tmp))
    {
        hB2_l = PyDouble_to_Double(hB2_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB2_r_tmp))
    {
        hB2_r = PyDouble_to_Double(hB2_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hPSI_l_tmp))
    {
        hPSI_l = PyDouble_to_Double(hPSI_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hPSI_r_tmp))
    {
        hPSI_r = PyDouble_to_Double(hPSI_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB1c_l_tmp))
    {
        hB1c_l = PyDouble_to_Double(hB1c_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB1c_r_tmp))
    {
        hB1c_r = PyDouble_to_Double(hB1c_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB2c_l_tmp))
    {
        hB2c_l = PyDouble_to_Double(hB2c_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(hB2c_r_tmp))
    {
        hB2c_r = PyDouble_to_Double(hB2c_r_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(Z_l_tmp))
    {
        Z_l = PyDouble_to_Double(Z_l_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(Z_r_tmp))
    {
        Z_r = PyDouble_to_Double(Z_r_tmp);
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
    if (PyIs_NativeFloat(mesure_tmp))
    {
        mesure = PyDouble_to_Double(mesure_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (PyIs_NativeFloat(grav_tmp))
    {
        grav = PyDouble_to_Double(grav_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(flux_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        flux = pyarray_to_ndarray(flux_tmp);
    }
    if (PyIs_NativeFloat(cpsi_tmp))
    {
        cpsi = PyDouble_to_Double(cpsi_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    bind_c_srnh_scheme_mhd(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, hPSI_l, hPSI_r, hB1c_l, hB1c_r, hB2c_l, hB2c_r, Z_l, Z_r, nd_ndim(&normal, 0), nd_data(&normal), mesure, grav, nd_ndim(&flux, 0), nd_data(&flux), cpsi);
    result = Py_BuildValue("");
    free_pointer(normal);
    free_pointer(flux);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *explicitscheme_convective_SWMHD_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray rez_h = {.shape = NULL};
    t_ndarray rez_hu = {.shape = NULL};
    t_ndarray rez_hv = {.shape = NULL};
    t_ndarray rez_hB1 = {.shape = NULL};
    t_ndarray rez_hB2 = {.shape = NULL};
    t_ndarray rez_PSI = {.shape = NULL};
    t_ndarray rez_Z = {.shape = NULL};
    t_ndarray h_c = {.shape = NULL};
    t_ndarray hu_c = {.shape = NULL};
    t_ndarray hv_c = {.shape = NULL};
    t_ndarray hB1_c = {.shape = NULL};
    t_ndarray hB2_c = {.shape = NULL};
    t_ndarray hPSIc = {.shape = NULL};
    t_ndarray Z_c = {.shape = NULL};
    t_ndarray h_ghost = {.shape = NULL};
    t_ndarray hu_ghost = {.shape = NULL};
    t_ndarray hv_ghost = {.shape = NULL};
    t_ndarray hB1_ghost = {.shape = NULL};
    t_ndarray hB2_ghost = {.shape = NULL};
    t_ndarray hPSIghost = {.shape = NULL};
    t_ndarray Z_ghost = {.shape = NULL};
    t_ndarray h_halo = {.shape = NULL};
    t_ndarray hu_halo = {.shape = NULL};
    t_ndarray hv_halo = {.shape = NULL};
    t_ndarray hB1_halo = {.shape = NULL};
    t_ndarray hB2_halo = {.shape = NULL};
    t_ndarray hPSIhalo = {.shape = NULL};
    t_ndarray Z_halo = {.shape = NULL};
    t_ndarray h_x = {.shape = NULL};
    t_ndarray h_y = {.shape = NULL};
    t_ndarray hx_halo = {.shape = NULL};
    t_ndarray hy_halo = {.shape = NULL};
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
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray boundaryfaces = {.shape = NULL};
    t_ndarray periodicboundaryfaces = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t order;
    double cpsi;
    t_ndarray hB1_cst = {.shape = NULL};
    t_ndarray hB2_cst = {.shape = NULL};
    PyArrayObject *rez_h_tmp;
    PyArrayObject *rez_hu_tmp;
    PyArrayObject *rez_hv_tmp;
    PyArrayObject *rez_hB1_tmp;
    PyArrayObject *rez_hB2_tmp;
    PyArrayObject *rez_PSI_tmp;
    PyArrayObject *rez_Z_tmp;
    PyArrayObject *h_c_tmp;
    PyArrayObject *hu_c_tmp;
    PyArrayObject *hv_c_tmp;
    PyArrayObject *hB1_c_tmp;
    PyArrayObject *hB2_c_tmp;
    PyArrayObject *hPSIc_tmp;
    PyArrayObject *Z_c_tmp;
    PyArrayObject *h_ghost_tmp;
    PyArrayObject *hu_ghost_tmp;
    PyArrayObject *hv_ghost_tmp;
    PyArrayObject *hB1_ghost_tmp;
    PyArrayObject *hB2_ghost_tmp;
    PyArrayObject *hPSIghost_tmp;
    PyArrayObject *Z_ghost_tmp;
    PyArrayObject *h_halo_tmp;
    PyArrayObject *hu_halo_tmp;
    PyArrayObject *hv_halo_tmp;
    PyArrayObject *hB1_halo_tmp;
    PyArrayObject *hB2_halo_tmp;
    PyArrayObject *hPSIhalo_tmp;
    PyArrayObject *Z_halo_tmp;
    PyArrayObject *h_x_tmp;
    PyArrayObject *h_y_tmp;
    PyArrayObject *hx_halo_tmp;
    PyArrayObject *hy_halo_tmp;
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
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *boundaryfaces_tmp;
    PyArrayObject *periodicboundaryfaces_tmp;
    PyArrayObject *shift_tmp;
    PyObject *order_tmp;
    PyObject *cpsi_tmp;
    PyArrayObject *hB1_cst_tmp;
    PyArrayObject *hB2_cst_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "rez_h",
        "rez_hu",
        "rez_hv",
        "rez_hB1",
        "rez_hB2",
        "rez_PSI",
        "rez_Z",
        "h_c",
        "hu_c",
        "hv_c",
        "hB1_c",
        "hB2_c",
        "hPSIc",
        "Z_c",
        "h_ghost",
        "hu_ghost",
        "hv_ghost",
        "hB1_ghost",
        "hB2_ghost",
        "hPSIghost",
        "Z_ghost",
        "h_halo",
        "hu_halo",
        "hv_halo",
        "hB1_halo",
        "hB2_halo",
        "hPSIhalo",
        "Z_halo",
        "h_x",
        "h_y",
        "hx_halo",
        "hy_halo",
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
        "innerfaces",
        "halofaces",
        "boundaryfaces",
        "periodicboundaryfaces",
        "shift",
        "order",
        "cpsi",
        "hB1_cst",
        "hB2_cst",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!OOO!O!", kwlist, &PyArray_Type, &rez_h_tmp, &PyArray_Type, &rez_hu_tmp, &PyArray_Type, &rez_hv_tmp, &PyArray_Type, &rez_hB1_tmp, &PyArray_Type, &rez_hB2_tmp, &PyArray_Type, &rez_PSI_tmp, &PyArray_Type, &rez_Z_tmp, &PyArray_Type, &h_c_tmp, &PyArray_Type, &hu_c_tmp, &PyArray_Type, &hv_c_tmp, &PyArray_Type, &hB1_c_tmp, &PyArray_Type, &hB2_c_tmp, &PyArray_Type, &hPSIc_tmp, &PyArray_Type, &Z_c_tmp, &PyArray_Type, &h_ghost_tmp, &PyArray_Type, &hu_ghost_tmp, &PyArray_Type, &hv_ghost_tmp, &PyArray_Type, &hB1_ghost_tmp, &PyArray_Type, &hB2_ghost_tmp, &PyArray_Type, &hPSIghost_tmp, &PyArray_Type, &Z_ghost_tmp, &PyArray_Type, &h_halo_tmp, &PyArray_Type, &hu_halo_tmp, &PyArray_Type, &hv_halo_tmp, &PyArray_Type, &hB1_halo_tmp, &PyArray_Type, &hB2_halo_tmp, &PyArray_Type, &hPSIhalo_tmp, &PyArray_Type, &Z_halo_tmp, &PyArray_Type, &h_x_tmp, &PyArray_Type, &h_y_tmp, &PyArray_Type, &hx_halo_tmp, &PyArray_Type, &hy_halo_tmp, &PyArray_Type, &psi_tmp, &PyArray_Type, &psi_halo_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerf_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &centerg_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &mesuref_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &boundaryfaces_tmp, &PyArray_Type, &periodicboundaryfaces_tmp, &PyArray_Type, &shift_tmp, &order_tmp, &cpsi_tmp, &PyArray_Type, &hB1_cst_tmp, &PyArray_Type, &hB2_cst_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(rez_h_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_h = pyarray_to_ndarray(rez_h_tmp);
    }
    if (!pyarray_check(rez_hu_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hu = pyarray_to_ndarray(rez_hu_tmp);
    }
    if (!pyarray_check(rez_hv_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hv = pyarray_to_ndarray(rez_hv_tmp);
    }
    if (!pyarray_check(rez_hB1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hB1 = pyarray_to_ndarray(rez_hB1_tmp);
    }
    if (!pyarray_check(rez_hB2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_hB2 = pyarray_to_ndarray(rez_hB2_tmp);
    }
    if (!pyarray_check(rez_PSI_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_PSI = pyarray_to_ndarray(rez_PSI_tmp);
    }
    if (!pyarray_check(rez_Z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_Z = pyarray_to_ndarray(rez_Z_tmp);
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(hu_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_c = pyarray_to_ndarray(hu_c_tmp);
    }
    if (!pyarray_check(hv_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_c = pyarray_to_ndarray(hv_c_tmp);
    }
    if (!pyarray_check(hB1_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_c = pyarray_to_ndarray(hB1_c_tmp);
    }
    if (!pyarray_check(hB2_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_c = pyarray_to_ndarray(hB2_c_tmp);
    }
    if (!pyarray_check(hPSIc_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hPSIc = pyarray_to_ndarray(hPSIc_tmp);
    }
    if (!pyarray_check(Z_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_c = pyarray_to_ndarray(Z_c_tmp);
    }
    if (!pyarray_check(h_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_ghost = pyarray_to_ndarray(h_ghost_tmp);
    }
    if (!pyarray_check(hu_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_ghost = pyarray_to_ndarray(hu_ghost_tmp);
    }
    if (!pyarray_check(hv_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_ghost = pyarray_to_ndarray(hv_ghost_tmp);
    }
    if (!pyarray_check(hB1_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_ghost = pyarray_to_ndarray(hB1_ghost_tmp);
    }
    if (!pyarray_check(hB2_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_ghost = pyarray_to_ndarray(hB2_ghost_tmp);
    }
    if (!pyarray_check(hPSIghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hPSIghost = pyarray_to_ndarray(hPSIghost_tmp);
    }
    if (!pyarray_check(Z_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_ghost = pyarray_to_ndarray(Z_ghost_tmp);
    }
    if (!pyarray_check(h_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_halo = pyarray_to_ndarray(h_halo_tmp);
    }
    if (!pyarray_check(hu_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hu_halo = pyarray_to_ndarray(hu_halo_tmp);
    }
    if (!pyarray_check(hv_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hv_halo = pyarray_to_ndarray(hv_halo_tmp);
    }
    if (!pyarray_check(hB1_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_halo = pyarray_to_ndarray(hB1_halo_tmp);
    }
    if (!pyarray_check(hB2_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_halo = pyarray_to_ndarray(hB2_halo_tmp);
    }
    if (!pyarray_check(hPSIhalo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hPSIhalo = pyarray_to_ndarray(hPSIhalo_tmp);
    }
    if (!pyarray_check(Z_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Z_halo = pyarray_to_ndarray(Z_halo_tmp);
    }
    if (!pyarray_check(h_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_x = pyarray_to_ndarray(h_x_tmp);
    }
    if (!pyarray_check(h_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_y = pyarray_to_ndarray(h_y_tmp);
    }
    if (!pyarray_check(hx_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hx_halo = pyarray_to_ndarray(hx_halo_tmp);
    }
    if (!pyarray_check(hy_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hy_halo = pyarray_to_ndarray(hy_halo_tmp);
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
    if (PyIs_NativeFloat(cpsi_tmp))
    {
        cpsi = PyDouble_to_Double(cpsi_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    if (!pyarray_check(hB1_cst_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB1_cst = pyarray_to_ndarray(hB1_cst_tmp);
    }
    if (!pyarray_check(hB2_cst_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        hB2_cst = pyarray_to_ndarray(hB2_cst_tmp);
    }
    bind_c_explicitscheme_convective_swmhd(nd_ndim(&rez_h, 0), nd_data(&rez_h), nd_ndim(&rez_hu, 0), nd_data(&rez_hu), nd_ndim(&rez_hv, 0), nd_data(&rez_hv), nd_ndim(&rez_hB1, 0), nd_data(&rez_hB1), nd_ndim(&rez_hB2, 0), nd_data(&rez_hB2), nd_ndim(&rez_PSI, 0), nd_data(&rez_PSI), nd_ndim(&rez_Z, 0), nd_data(&rez_Z), nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&hu_c, 0), nd_data(&hu_c), nd_ndim(&hv_c, 0), nd_data(&hv_c), nd_ndim(&hB1_c, 0), nd_data(&hB1_c), nd_ndim(&hB2_c, 0), nd_data(&hB2_c), nd_ndim(&hPSIc, 0), nd_data(&hPSIc), nd_ndim(&Z_c, 0), nd_data(&Z_c), nd_ndim(&h_ghost, 0), nd_data(&h_ghost), nd_ndim(&hu_ghost, 0), nd_data(&hu_ghost), nd_ndim(&hv_ghost, 0), nd_data(&hv_ghost), nd_ndim(&hB1_ghost, 0), nd_data(&hB1_ghost), nd_ndim(&hB2_ghost, 0), nd_data(&hB2_ghost), nd_ndim(&hPSIghost, 0), nd_data(&hPSIghost), nd_ndim(&Z_ghost, 0), nd_data(&Z_ghost), nd_ndim(&h_halo, 0), nd_data(&h_halo), nd_ndim(&hu_halo, 0), nd_data(&hu_halo), nd_ndim(&hv_halo, 0), nd_data(&hv_halo), nd_ndim(&hB1_halo, 0), nd_data(&hB1_halo), nd_ndim(&hB2_halo, 0), nd_data(&hB2_halo), nd_ndim(&hPSIhalo, 0), nd_data(&hPSIhalo), nd_ndim(&Z_halo, 0), nd_data(&Z_halo), nd_ndim(&h_x, 0), nd_data(&h_x), nd_ndim(&h_y, 0), nd_data(&h_y), nd_ndim(&hx_halo, 0), nd_data(&hx_halo), nd_ndim(&hy_halo, 0), nd_data(&hy_halo), nd_ndim(&psi, 0), nd_data(&psi), nd_ndim(&psi_halo, 0), nd_data(&psi_halo), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&centerg, 0), nd_ndim(&centerg, 1), nd_data(&centerg), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&mesuref, 0), nd_data(&mesuref), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&boundaryfaces, 0), nd_data(&boundaryfaces), nd_ndim(&periodicboundaryfaces, 0), nd_data(&periodicboundaryfaces), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), order, cpsi, nd_ndim(&hB1_cst, 0), nd_data(&hB1_cst), nd_ndim(&hB2_cst, 0), nd_data(&hB2_cst));
    result = Py_BuildValue("");
    free_pointer(rez_h);
    free_pointer(rez_hu);
    free_pointer(rez_hv);
    free_pointer(rez_hB1);
    free_pointer(rez_hB2);
    free_pointer(rez_PSI);
    free_pointer(rez_Z);
    free_pointer(h_c);
    free_pointer(hu_c);
    free_pointer(hv_c);
    free_pointer(hB1_c);
    free_pointer(hB2_c);
    free_pointer(hPSIc);
    free_pointer(Z_c);
    free_pointer(h_ghost);
    free_pointer(hu_ghost);
    free_pointer(hv_ghost);
    free_pointer(hB1_ghost);
    free_pointer(hB2_ghost);
    free_pointer(hPSIghost);
    free_pointer(Z_ghost);
    free_pointer(h_halo);
    free_pointer(hu_halo);
    free_pointer(hv_halo);
    free_pointer(hB1_halo);
    free_pointer(hB2_halo);
    free_pointer(hPSIhalo);
    free_pointer(Z_halo);
    free_pointer(h_x);
    free_pointer(h_y);
    free_pointer(hx_halo);
    free_pointer(hy_halo);
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
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(boundaryfaces);
    free_pointer(periodicboundaryfaces);
    free_pointer(shift);
    free_pointer(hB1_cst);
    free_pointer(hB2_cst);
    return result;
}
/*........................................*/

static int exec_func(PyObject* m)
{
    return 0;
}

/*........................................*/

static PyMethodDef tools_SWMHD_methods[] = {
    {
        "Total_Energy",
        (PyCFunction)Total_Energy_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "explicitscheme_stabilizater",
        (PyCFunction)explicitscheme_stabilizater_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "time_step_SWMHD",
        (PyCFunction)time_step_SWMHD_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "update_SWMHD",
        (PyCFunction)update_SWMHD_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "cpsi_global",
        (PyCFunction)cpsi_global_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "initialisation_SWMHD",
        (PyCFunction)initialisation_SWMHD_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "term_source_srnh_SWMHD",
        (PyCFunction)term_source_srnh_SWMHD_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "LF_scheme_MHD",
        (PyCFunction)LF_scheme_MHD_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "srnh_scheme_MHD",
        (PyCFunction)srnh_scheme_MHD_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "explicitscheme_convective_SWMHD",
        (PyCFunction)explicitscheme_convective_SWMHD_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static PyModuleDef_Slot tools_SWMHD_slots[] = {
    {Py_mod_exec, exec_func},
    {0, NULL},
};

/*........................................*/

static struct PyModuleDef tools_SWMHD_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "tools_SWMHD",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    tools_SWMHD_methods,
    tools_SWMHD_slots
};

/*........................................*/

PyMODINIT_FUNC PyInit_tools_SWMHD(void)
{
    import_array();
    return PyModuleDef_Init(&tools_SWMHD_module);
}
