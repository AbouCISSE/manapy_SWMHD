#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include <stdlib.h>
#include <stdint.h>
#include "ndarrays.h"
#include "cwrapper_ndarrays.h"


void bind_c_compute_k(int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_namef, int64_t *namef, int64_t n0_ghostcenterf, int64_t n1_ghostcenterf, double *ghostcenterf, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_K, int64_t n1_K, double *K, int64_t dim);
void bind_c_create_info_2dfaces(int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_namen, int64_t *namen, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t nbfaces, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_mesuref, double *mesuref, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_namef, int64_t *namef);
void bind_c_create_info_3dfaces(int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_namen, int64_t *namen, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t nbfaces, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_mesuref, double *mesuref, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_namef, int64_t *namef);
void bind_c_compute_2dcentervolumeofcell(int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t nbelements, int64_t n0_center, int64_t n1_center, double *center, int64_t n0_volume, double *volume);
void bind_c_compute_3dcentervolumeofcell(int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t nbelements, int64_t n0_center, int64_t n1_center, double *center, int64_t n0_volume, double *volume);
void bind_c_create_cellsofface(int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t nbelements, int64_t nbfaces, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t dim);
void bind_c_create_2dfaces(int64_t n0_nodeidc, int64_t n1_nodeidc, int64_t *nodeidc, int64_t nbelements, int64_t n0_faces, int64_t n1_faces, int64_t *faces, int64_t n0_cellf, int64_t n1_cellf, int64_t *cellf);
void bind_c_create_cell_faceid(int64_t nbelements, int64_t n0_oldTonewIndex, int64_t *oldTonewIndex, int64_t n0_cellf, int64_t n1_cellf, int64_t *cellf, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t dim);
void bind_c_create_3dfaces(int64_t n0_nodeidc, int64_t n1_nodeidc, int64_t *nodeidc, int64_t nbelements, int64_t n0_faces, int64_t n1_faces, int64_t *faces, int64_t n0_cellf, int64_t n1_cellf, int64_t *cellf);
void bind_c_create_normalfacesofcell(int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t n0_normal, int64_t n1_normal, double *normal, int64_t nbelements, int64_t n0_nf, int64_t n1_nf, int64_t n2_nf, double *nf, int64_t dim);
void bind_c_face_gradient_info_2d(int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_namef, int64_t *namef, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_halofid, int64_t *halofid, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_airDiamond, double *airDiamond, int64_t n0_param1, double *param1, int64_t n0_param2, double *param2, int64_t n0_param3, double *param3, int64_t n0_param4, double *param4, int64_t n0_f_1, int64_t n1_f_1, double *f_1, int64_t n0_f_2, int64_t n1_f_2, double *f_2, int64_t n0_f_3, int64_t n1_f_3, double *f_3, int64_t n0_f_4, int64_t n1_f_4, double *f_4, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t dim);
void bind_c_variables(int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_cellidn, int64_t n1_cellidn, int64_t *cellidn, int64_t n0_haloidn, int64_t n1_haloidn, int64_t *haloidn, int64_t n0_periodicn, int64_t n1_periodicn, int64_t *periodicn, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_namen, int64_t *namen, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocentergn, int64_t n1_halocentergn, int64_t n2_halocentergn, double *halocentergn, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t nbproc, int64_t n0_R_x, double *R_x, int64_t n0_R_y, double *R_y, int64_t n0_lambda_x, double *lambda_x, int64_t n0_lambda_y, double *lambda_y, int64_t n0_number, int64_t *number, int64_t n0_shift, int64_t n1_shift, double *shift);
void bind_c_variables_3d(int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_cellidn, int64_t n1_cellidn, int64_t *cellidn, int64_t n0_haloidn, int64_t n1_haloidn, int64_t *haloidn, int64_t n0_periodicn, int64_t n1_periodicn, int64_t *periodicn, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_namen, int64_t *namen, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocenterg, int64_t n1_halocenterg, int64_t n2_halocenterg, double *halocenterg, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t nbproc, int64_t n0_R_x, double *R_x, int64_t n0_R_y, double *R_y, int64_t n0_R_z, double *R_z, int64_t n0_lambda_x, double *lambda_x, int64_t n0_lambda_y, double *lambda_y, int64_t n0_lambda_z, double *lambda_z, int64_t n0_number, int64_t *number, int64_t n0_shift, int64_t n1_shift, double *shift);

/*........................................*/


/*........................................*/

/*........................................*/
PyObject *compute_K_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellid = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray ghostcenterf = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray K = {.shape = NULL};
    int64_t dim;
    PyArrayObject *cellid_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *ghostcenterf_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *K_tmp;
    PyObject *dim_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellid",
        "namef",
        "ghostcenterf",
        "centerc",
        "K",
        "dim",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O", kwlist, &PyArray_Type, &cellid_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &ghostcenterf_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &K_tmp, &dim_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(ghostcenterf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        ghostcenterf = pyarray_to_ndarray(ghostcenterf_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(K_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        K = pyarray_to_ndarray(K_tmp);
    }
    if (PyIs_NativeInt(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_compute_k(nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&ghostcenterf, 0), nd_ndim(&ghostcenterf, 1), nd_data(&ghostcenterf), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&K, 0), nd_ndim(&K, 1), nd_data(&K), dim);
    result = Py_BuildValue("");
    free_pointer(cellid);
    free_pointer(namef);
    free_pointer(ghostcenterf);
    free_pointer(centerc);
    free_pointer(K);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_info_2dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellid = {.shape = NULL};
    t_ndarray nodeid = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray vertex = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    int64_t nbfaces;
    t_ndarray normalf = {.shape = NULL};
    t_ndarray mesuref = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    PyArrayObject *cellid_tmp;
    PyArrayObject *nodeid_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *vertex_tmp;
    PyArrayObject *centerc_tmp;
    PyObject *nbfaces_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *mesuref_tmp;
    PyArrayObject *centerf_tmp;
    PyArrayObject *namef_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellid",
        "nodeid",
        "namen",
        "vertex",
        "centerc",
        "nbfaces",
        "normalf",
        "mesuref",
        "centerf",
        "namef",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!OO!O!O!O!", kwlist, &PyArray_Type, &cellid_tmp, &PyArray_Type, &nodeid_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &vertex_tmp, &PyArray_Type, &centerc_tmp, &nbfaces_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &mesuref_tmp, &PyArray_Type, &centerf_tmp, &PyArray_Type, &namef_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = pyarray_to_ndarray(nodeid_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = pyarray_to_ndarray(vertex_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (PyIs_NativeInt(nbfaces_tmp))
    {
        nbfaces = PyInt64_to_Int64(nbfaces_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(mesuref_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesuref = pyarray_to_ndarray(mesuref_tmp);
    }
    if (!pyarray_check(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = pyarray_to_ndarray(centerf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    bind_c_create_info_2dfaces(nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&nodeid, 0), nd_ndim(&nodeid, 1), nd_data(&nodeid), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&vertex, 0), nd_ndim(&vertex, 1), nd_data(&vertex), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nbfaces, nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&mesuref, 0), nd_data(&mesuref), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf), nd_ndim(&namef, 0), nd_data(&namef));
    result = Py_BuildValue("");
    free_pointer(cellid);
    free_pointer(nodeid);
    free_pointer(namen);
    free_pointer(vertex);
    free_pointer(centerc);
    free_pointer(normalf);
    free_pointer(mesuref);
    free_pointer(centerf);
    free_pointer(namef);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_info_3dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellid = {.shape = NULL};
    t_ndarray nodeid = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray vertex = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    int64_t nbfaces;
    t_ndarray normalf = {.shape = NULL};
    t_ndarray mesuref = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    PyArrayObject *cellid_tmp;
    PyArrayObject *nodeid_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *vertex_tmp;
    PyArrayObject *centerc_tmp;
    PyObject *nbfaces_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *mesuref_tmp;
    PyArrayObject *centerf_tmp;
    PyArrayObject *namef_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellid",
        "nodeid",
        "namen",
        "vertex",
        "centerc",
        "nbfaces",
        "normalf",
        "mesuref",
        "centerf",
        "namef",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!OO!O!O!O!", kwlist, &PyArray_Type, &cellid_tmp, &PyArray_Type, &nodeid_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &vertex_tmp, &PyArray_Type, &centerc_tmp, &nbfaces_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &mesuref_tmp, &PyArray_Type, &centerf_tmp, &PyArray_Type, &namef_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = pyarray_to_ndarray(nodeid_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = pyarray_to_ndarray(vertex_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (PyIs_NativeInt(nbfaces_tmp))
    {
        nbfaces = PyInt64_to_Int64(nbfaces_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(mesuref_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesuref = pyarray_to_ndarray(mesuref_tmp);
    }
    if (!pyarray_check(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = pyarray_to_ndarray(centerf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    bind_c_create_info_3dfaces(nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&nodeid, 0), nd_ndim(&nodeid, 1), nd_data(&nodeid), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&vertex, 0), nd_ndim(&vertex, 1), nd_data(&vertex), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nbfaces, nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&mesuref, 0), nd_data(&mesuref), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf), nd_ndim(&namef, 0), nd_data(&namef));
    result = Py_BuildValue("");
    free_pointer(cellid);
    free_pointer(nodeid);
    free_pointer(namen);
    free_pointer(vertex);
    free_pointer(centerc);
    free_pointer(normalf);
    free_pointer(mesuref);
    free_pointer(centerf);
    free_pointer(namef);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *Compute_2dcentervolumeOfCell_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray nodeid = {.shape = NULL};
    t_ndarray vertex = {.shape = NULL};
    int64_t nbelements;
    t_ndarray center = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    PyArrayObject *nodeid_tmp;
    PyArrayObject *vertex_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *center_tmp;
    PyArrayObject *volume_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "nodeid",
        "vertex",
        "nbelements",
        "center",
        "volume",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!OO!O!", kwlist, &PyArray_Type, &nodeid_tmp, &PyArray_Type, &vertex_tmp, &nbelements_tmp, &PyArray_Type, &center_tmp, &PyArray_Type, &volume_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = pyarray_to_ndarray(nodeid_tmp);
    }
    if (!pyarray_check(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = pyarray_to_ndarray(vertex_tmp);
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(center_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        center = pyarray_to_ndarray(center_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    bind_c_compute_2dcentervolumeofcell(nd_ndim(&nodeid, 0), nd_ndim(&nodeid, 1), nd_data(&nodeid), nd_ndim(&vertex, 0), nd_ndim(&vertex, 1), nd_data(&vertex), nbelements, nd_ndim(&center, 0), nd_ndim(&center, 1), nd_data(&center), nd_ndim(&volume, 0), nd_data(&volume));
    result = Py_BuildValue("");
    free_pointer(nodeid);
    free_pointer(vertex);
    free_pointer(center);
    free_pointer(volume);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *Compute_3dcentervolumeOfCell_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray nodeid = {.shape = NULL};
    t_ndarray vertex = {.shape = NULL};
    int64_t nbelements;
    t_ndarray center = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    PyArrayObject *nodeid_tmp;
    PyArrayObject *vertex_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *center_tmp;
    PyArrayObject *volume_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "nodeid",
        "vertex",
        "nbelements",
        "center",
        "volume",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!OO!O!", kwlist, &PyArray_Type, &nodeid_tmp, &PyArray_Type, &vertex_tmp, &nbelements_tmp, &PyArray_Type, &center_tmp, &PyArray_Type, &volume_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = pyarray_to_ndarray(nodeid_tmp);
    }
    if (!pyarray_check(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = pyarray_to_ndarray(vertex_tmp);
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(center_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        center = pyarray_to_ndarray(center_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    bind_c_compute_3dcentervolumeofcell(nd_ndim(&nodeid, 0), nd_ndim(&nodeid, 1), nd_data(&nodeid), nd_ndim(&vertex, 0), nd_ndim(&vertex, 1), nd_data(&vertex), nbelements, nd_ndim(&center, 0), nd_ndim(&center, 1), nd_data(&center), nd_ndim(&volume, 0), nd_data(&volume));
    result = Py_BuildValue("");
    free_pointer(nodeid);
    free_pointer(vertex);
    free_pointer(center);
    free_pointer(volume);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_cellsOfFace_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray faceid = {.shape = NULL};
    int64_t nbelements;
    int64_t nbfaces;
    t_ndarray cellid = {.shape = NULL};
    int64_t dim;
    PyArrayObject *faceid_tmp;
    PyObject *nbelements_tmp;
    PyObject *nbfaces_tmp;
    PyArrayObject *cellid_tmp;
    PyObject *dim_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "faceid",
        "nbelements",
        "nbfaces",
        "cellid",
        "dim",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!OOO!O", kwlist, &PyArray_Type, &faceid_tmp, &nbelements_tmp, &nbfaces_tmp, &PyArray_Type, &cellid_tmp, &dim_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = pyarray_to_ndarray(faceid_tmp);
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (PyIs_NativeInt(nbfaces_tmp))
    {
        nbfaces = PyInt64_to_Int64(nbfaces_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (PyIs_NativeInt(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_create_cellsofface(nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid), nbelements, nbfaces, nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), dim);
    result = Py_BuildValue("");
    free_pointer(faceid);
    free_pointer(cellid);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_2dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray nodeidc = {.shape = NULL};
    int64_t nbelements;
    t_ndarray faces = {.shape = NULL};
    t_ndarray cellf = {.shape = NULL};
    PyArrayObject *nodeidc_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *faces_tmp;
    PyArrayObject *cellf_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "nodeidc",
        "nbelements",
        "faces",
        "cellf",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!OO!O!", kwlist, &PyArray_Type, &nodeidc_tmp, &nbelements_tmp, &PyArray_Type, &faces_tmp, &PyArray_Type, &cellf_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(nodeidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidc = pyarray_to_ndarray(nodeidc_tmp);
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(faces_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faces = pyarray_to_ndarray(faces_tmp);
    }
    if (!pyarray_check(cellf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellf = pyarray_to_ndarray(cellf_tmp);
    }
    bind_c_create_2dfaces(nd_ndim(&nodeidc, 0), nd_ndim(&nodeidc, 1), nd_data(&nodeidc), nbelements, nd_ndim(&faces, 0), nd_ndim(&faces, 1), nd_data(&faces), nd_ndim(&cellf, 0), nd_ndim(&cellf, 1), nd_data(&cellf));
    result = Py_BuildValue("");
    free_pointer(nodeidc);
    free_pointer(faces);
    free_pointer(cellf);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_cell_faceid_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t nbelements;
    t_ndarray oldTonewIndex = {.shape = NULL};
    t_ndarray cellf = {.shape = NULL};
    t_ndarray faceid = {.shape = NULL};
    int64_t dim;
    PyObject *nbelements_tmp;
    PyArrayObject *oldTonewIndex_tmp;
    PyArrayObject *cellf_tmp;
    PyArrayObject *faceid_tmp;
    PyObject *dim_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "nbelements",
        "oldTonewIndex",
        "cellf",
        "faceid",
        "dim",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO!O!O!O", kwlist, &nbelements_tmp, &PyArray_Type, &oldTonewIndex_tmp, &PyArray_Type, &cellf_tmp, &PyArray_Type, &faceid_tmp, &dim_tmp))
    {
        return NULL;
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(oldTonewIndex_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldTonewIndex = pyarray_to_ndarray(oldTonewIndex_tmp);
    }
    if (!pyarray_check(cellf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellf = pyarray_to_ndarray(cellf_tmp);
    }
    if (!pyarray_check(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = pyarray_to_ndarray(faceid_tmp);
    }
    if (PyIs_NativeInt(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_create_cell_faceid(nbelements, nd_ndim(&oldTonewIndex, 0), nd_data(&oldTonewIndex), nd_ndim(&cellf, 0), nd_ndim(&cellf, 1), nd_data(&cellf), nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid), dim);
    result = Py_BuildValue("");
    free_pointer(oldTonewIndex);
    free_pointer(cellf);
    free_pointer(faceid);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_3dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray nodeidc = {.shape = NULL};
    int64_t nbelements;
    t_ndarray faces = {.shape = NULL};
    t_ndarray cellf = {.shape = NULL};
    PyArrayObject *nodeidc_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *faces_tmp;
    PyArrayObject *cellf_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "nodeidc",
        "nbelements",
        "faces",
        "cellf",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!OO!O!", kwlist, &PyArray_Type, &nodeidc_tmp, &nbelements_tmp, &PyArray_Type, &faces_tmp, &PyArray_Type, &cellf_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(nodeidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidc = pyarray_to_ndarray(nodeidc_tmp);
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(faces_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faces = pyarray_to_ndarray(faces_tmp);
    }
    if (!pyarray_check(cellf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellf = pyarray_to_ndarray(cellf_tmp);
    }
    bind_c_create_3dfaces(nd_ndim(&nodeidc, 0), nd_ndim(&nodeidc, 1), nd_data(&nodeidc), nbelements, nd_ndim(&faces, 0), nd_ndim(&faces, 1), nd_data(&faces), nd_ndim(&cellf, 0), nd_ndim(&cellf, 1), nd_data(&cellf));
    result = Py_BuildValue("");
    free_pointer(nodeidc);
    free_pointer(faces);
    free_pointer(cellf);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_NormalFacesOfCell_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    t_ndarray faceid = {.shape = NULL};
    t_ndarray normal = {.shape = NULL};
    int64_t nbelements;
    t_ndarray nf = {.shape = NULL};
    int64_t dim;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerf_tmp;
    PyArrayObject *faceid_tmp;
    PyArrayObject *normal_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *nf_tmp;
    PyObject *dim_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "centerc",
        "centerf",
        "faceid",
        "normal",
        "nbelements",
        "nf",
        "dim",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!OO!O", kwlist, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerf_tmp, &PyArray_Type, &faceid_tmp, &PyArray_Type, &normal_tmp, &nbelements_tmp, &PyArray_Type, &nf_tmp, &dim_tmp))
    {
        return NULL;
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
    if (!pyarray_check(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = pyarray_to_ndarray(faceid_tmp);
    }
    if (!pyarray_check(normal_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normal = pyarray_to_ndarray(normal_tmp);
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(nf_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nf = pyarray_to_ndarray(nf_tmp);
    }
    if (PyIs_NativeInt(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_create_normalfacesofcell(nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf), nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid), nd_ndim(&normal, 0), nd_ndim(&normal, 1), nd_data(&normal), nbelements, nd_ndim(&nf, 0), nd_ndim(&nf, 1), nd_ndim(&nf, 2), nd_data(&nf), dim);
    result = Py_BuildValue("");
    free_pointer(centerc);
    free_pointer(centerf);
    free_pointer(faceid);
    free_pointer(normal);
    free_pointer(nf);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *face_gradient_info_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray param1 = {.shape = NULL};
    t_ndarray param2 = {.shape = NULL};
    t_ndarray param3 = {.shape = NULL};
    t_ndarray param4 = {.shape = NULL};
    t_ndarray f_1 = {.shape = NULL};
    t_ndarray f_2 = {.shape = NULL};
    t_ndarray f_3 = {.shape = NULL};
    t_ndarray f_4 = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t dim;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *param1_tmp;
    PyArrayObject *param2_tmp;
    PyArrayObject *param3_tmp;
    PyArrayObject *param4_tmp;
    PyArrayObject *f_1_tmp;
    PyArrayObject *f_2_tmp;
    PyArrayObject *f_3_tmp;
    PyArrayObject *f_4_tmp;
    PyArrayObject *shift_tmp;
    PyObject *dim_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellidf",
        "nodeidf",
        "centergf",
        "namef",
        "normalf",
        "centerc",
        "centerh",
        "halofid",
        "vertexn",
        "airDiamond",
        "param1",
        "param2",
        "param3",
        "param4",
        "f_1",
        "f_2",
        "f_3",
        "f_4",
        "shift",
        "dim",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O", kwlist, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &param1_tmp, &PyArray_Type, &param2_tmp, &PyArray_Type, &param3_tmp, &PyArray_Type, &param4_tmp, &PyArray_Type, &f_1_tmp, &PyArray_Type, &f_2_tmp, &PyArray_Type, &f_3_tmp, &PyArray_Type, &f_4_tmp, &PyArray_Type, &shift_tmp, &dim_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(param1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param1 = pyarray_to_ndarray(param1_tmp);
    }
    if (!pyarray_check(param2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param2 = pyarray_to_ndarray(param2_tmp);
    }
    if (!pyarray_check(param3_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param3 = pyarray_to_ndarray(param3_tmp);
    }
    if (!pyarray_check(param4_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param4 = pyarray_to_ndarray(param4_tmp);
    }
    if (!pyarray_check(f_1_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_1 = pyarray_to_ndarray(f_1_tmp);
    }
    if (!pyarray_check(f_2_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_2 = pyarray_to_ndarray(f_2_tmp);
    }
    if (!pyarray_check(f_3_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_3 = pyarray_to_ndarray(f_3_tmp);
    }
    if (!pyarray_check(f_4_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_4 = pyarray_to_ndarray(f_4_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_face_gradient_info_2d(nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&param1, 0), nd_data(&param1), nd_ndim(&param2, 0), nd_data(&param2), nd_ndim(&param3, 0), nd_data(&param3), nd_ndim(&param4, 0), nd_data(&param4), nd_ndim(&f_1, 0), nd_ndim(&f_1, 1), nd_data(&f_1), nd_ndim(&f_2, 0), nd_ndim(&f_2, 1), nd_data(&f_2), nd_ndim(&f_3, 0), nd_ndim(&f_3, 1), nd_data(&f_3), nd_ndim(&f_4, 0), nd_ndim(&f_4, 1), nd_data(&f_4), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), dim);
    result = Py_BuildValue("");
    free_pointer(cellidf);
    free_pointer(nodeidf);
    free_pointer(centergf);
    free_pointer(namef);
    free_pointer(normalf);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(halofid);
    free_pointer(vertexn);
    free_pointer(airDiamond);
    free_pointer(param1);
    free_pointer(param2);
    free_pointer(param3);
    free_pointer(param4);
    free_pointer(f_1);
    free_pointer(f_2);
    free_pointer(f_3);
    free_pointer(f_4);
    free_pointer(shift);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *variables_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray centerc = {.shape = NULL};
    t_ndarray cellidn = {.shape = NULL};
    t_ndarray haloidn = {.shape = NULL};
    t_ndarray periodicn = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocentergn = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    int64_t nbproc;
    t_ndarray R_x = {.shape = NULL};
    t_ndarray R_y = {.shape = NULL};
    t_ndarray lambda_x = {.shape = NULL};
    t_ndarray lambda_y = {.shape = NULL};
    t_ndarray number = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    PyArrayObject *centerc_tmp;
    PyArrayObject *cellidn_tmp;
    PyArrayObject *haloidn_tmp;
    PyArrayObject *periodicn_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocentergn_tmp;
    PyArrayObject *centerh_tmp;
    PyObject *nbproc_tmp;
    PyArrayObject *R_x_tmp;
    PyArrayObject *R_y_tmp;
    PyArrayObject *lambda_x_tmp;
    PyArrayObject *lambda_y_tmp;
    PyArrayObject *number_tmp;
    PyArrayObject *shift_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "centerc",
        "cellidn",
        "haloidn",
        "periodicn",
        "vertexn",
        "namen",
        "centergn",
        "halocentergn",
        "centerh",
        "nbproc",
        "R_x",
        "R_y",
        "lambda_x",
        "lambda_y",
        "number",
        "shift",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!OO!O!O!O!O!O!", kwlist, &PyArray_Type, &centerc_tmp, &PyArray_Type, &cellidn_tmp, &PyArray_Type, &haloidn_tmp, &PyArray_Type, &periodicn_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocentergn_tmp, &PyArray_Type, &centerh_tmp, &nbproc_tmp, &PyArray_Type, &R_x_tmp, &PyArray_Type, &R_y_tmp, &PyArray_Type, &lambda_x_tmp, &PyArray_Type, &lambda_y_tmp, &PyArray_Type, &number_tmp, &PyArray_Type, &shift_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(cellidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidn = pyarray_to_ndarray(cellidn_tmp);
    }
    if (!pyarray_check(haloidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloidn = pyarray_to_ndarray(haloidn_tmp);
    }
    if (!pyarray_check(periodicn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicn = pyarray_to_ndarray(periodicn_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocentergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocentergn = pyarray_to_ndarray(halocentergn_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (PyIs_NativeInt(nbproc_tmp))
    {
        nbproc = PyInt64_to_Int64(nbproc_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(R_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_x = pyarray_to_ndarray(R_x_tmp);
    }
    if (!pyarray_check(R_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_y = pyarray_to_ndarray(R_y_tmp);
    }
    if (!pyarray_check(lambda_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_x = pyarray_to_ndarray(lambda_x_tmp);
    }
    if (!pyarray_check(lambda_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_y = pyarray_to_ndarray(lambda_y_tmp);
    }
    if (!pyarray_check(number_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        number = pyarray_to_ndarray(number_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    bind_c_variables(nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&cellidn, 0), nd_ndim(&cellidn, 1), nd_data(&cellidn), nd_ndim(&haloidn, 0), nd_ndim(&haloidn, 1), nd_data(&haloidn), nd_ndim(&periodicn, 0), nd_ndim(&periodicn, 1), nd_data(&periodicn), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocentergn, 0), nd_ndim(&halocentergn, 1), nd_ndim(&halocentergn, 2), nd_data(&halocentergn), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nbproc, nd_ndim(&R_x, 0), nd_data(&R_x), nd_ndim(&R_y, 0), nd_data(&R_y), nd_ndim(&lambda_x, 0), nd_data(&lambda_x), nd_ndim(&lambda_y, 0), nd_data(&lambda_y), nd_ndim(&number, 0), nd_data(&number), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift));
    result = Py_BuildValue("");
    free_pointer(centerc);
    free_pointer(cellidn);
    free_pointer(haloidn);
    free_pointer(periodicn);
    free_pointer(vertexn);
    free_pointer(namen);
    free_pointer(centergn);
    free_pointer(halocentergn);
    free_pointer(centerh);
    free_pointer(R_x);
    free_pointer(R_y);
    free_pointer(lambda_x);
    free_pointer(lambda_y);
    free_pointer(number);
    free_pointer(shift);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *variables_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray centerc = {.shape = NULL};
    t_ndarray cellidn = {.shape = NULL};
    t_ndarray haloidn = {.shape = NULL};
    t_ndarray periodicn = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocenterg = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    int64_t nbproc;
    t_ndarray R_x = {.shape = NULL};
    t_ndarray R_y = {.shape = NULL};
    t_ndarray R_z = {.shape = NULL};
    t_ndarray lambda_x = {.shape = NULL};
    t_ndarray lambda_y = {.shape = NULL};
    t_ndarray lambda_z = {.shape = NULL};
    t_ndarray number = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    PyArrayObject *centerc_tmp;
    PyArrayObject *cellidn_tmp;
    PyArrayObject *haloidn_tmp;
    PyArrayObject *periodicn_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocenterg_tmp;
    PyArrayObject *centerh_tmp;
    PyObject *nbproc_tmp;
    PyArrayObject *R_x_tmp;
    PyArrayObject *R_y_tmp;
    PyArrayObject *R_z_tmp;
    PyArrayObject *lambda_x_tmp;
    PyArrayObject *lambda_y_tmp;
    PyArrayObject *lambda_z_tmp;
    PyArrayObject *number_tmp;
    PyArrayObject *shift_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "centerc",
        "cellidn",
        "haloidn",
        "periodicn",
        "vertexn",
        "namen",
        "centergn",
        "halocenterg",
        "centerh",
        "nbproc",
        "R_x",
        "R_y",
        "R_z",
        "lambda_x",
        "lambda_y",
        "lambda_z",
        "number",
        "shift",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!OO!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &centerc_tmp, &PyArray_Type, &cellidn_tmp, &PyArray_Type, &haloidn_tmp, &PyArray_Type, &periodicn_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocenterg_tmp, &PyArray_Type, &centerh_tmp, &nbproc_tmp, &PyArray_Type, &R_x_tmp, &PyArray_Type, &R_y_tmp, &PyArray_Type, &R_z_tmp, &PyArray_Type, &lambda_x_tmp, &PyArray_Type, &lambda_y_tmp, &PyArray_Type, &lambda_z_tmp, &PyArray_Type, &number_tmp, &PyArray_Type, &shift_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(cellidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidn = pyarray_to_ndarray(cellidn_tmp);
    }
    if (!pyarray_check(haloidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloidn = pyarray_to_ndarray(haloidn_tmp);
    }
    if (!pyarray_check(periodicn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicn = pyarray_to_ndarray(periodicn_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocenterg_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocenterg = pyarray_to_ndarray(halocenterg_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (PyIs_NativeInt(nbproc_tmp))
    {
        nbproc = PyInt64_to_Int64(nbproc_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(R_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_x = pyarray_to_ndarray(R_x_tmp);
    }
    if (!pyarray_check(R_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_y = pyarray_to_ndarray(R_y_tmp);
    }
    if (!pyarray_check(R_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_z = pyarray_to_ndarray(R_z_tmp);
    }
    if (!pyarray_check(lambda_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_x = pyarray_to_ndarray(lambda_x_tmp);
    }
    if (!pyarray_check(lambda_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_y = pyarray_to_ndarray(lambda_y_tmp);
    }
    if (!pyarray_check(lambda_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_z = pyarray_to_ndarray(lambda_z_tmp);
    }
    if (!pyarray_check(number_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        number = pyarray_to_ndarray(number_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    bind_c_variables_3d(nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&cellidn, 0), nd_ndim(&cellidn, 1), nd_data(&cellidn), nd_ndim(&haloidn, 0), nd_ndim(&haloidn, 1), nd_data(&haloidn), nd_ndim(&periodicn, 0), nd_ndim(&periodicn, 1), nd_data(&periodicn), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocenterg, 0), nd_ndim(&halocenterg, 1), nd_ndim(&halocenterg, 2), nd_data(&halocenterg), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nbproc, nd_ndim(&R_x, 0), nd_data(&R_x), nd_ndim(&R_y, 0), nd_data(&R_y), nd_ndim(&R_z, 0), nd_data(&R_z), nd_ndim(&lambda_x, 0), nd_data(&lambda_x), nd_ndim(&lambda_y, 0), nd_data(&lambda_y), nd_ndim(&lambda_z, 0), nd_data(&lambda_z), nd_ndim(&number, 0), nd_data(&number), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift));
    result = Py_BuildValue("");
    free_pointer(centerc);
    free_pointer(cellidn);
    free_pointer(haloidn);
    free_pointer(periodicn);
    free_pointer(vertexn);
    free_pointer(namen);
    free_pointer(centergn);
    free_pointer(halocenterg);
    free_pointer(centerh);
    free_pointer(R_x);
    free_pointer(R_y);
    free_pointer(R_z);
    free_pointer(lambda_x);
    free_pointer(lambda_y);
    free_pointer(lambda_z);
    free_pointer(number);
    free_pointer(shift);
    return result;
}
/*........................................*/

static int exec_func(PyObject* m)
{
    return 0;
}

/*........................................*/

static PyMethodDef pyccel_ddm_methods[] = {
    {
        "compute_K",
        (PyCFunction)compute_K_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "create_info_2dfaces",
        (PyCFunction)create_info_2dfaces_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "create_info_3dfaces",
        (PyCFunction)create_info_3dfaces_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "Compute_2dcentervolumeOfCell",
        (PyCFunction)Compute_2dcentervolumeOfCell_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "Compute_3dcentervolumeOfCell",
        (PyCFunction)Compute_3dcentervolumeOfCell_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "create_cellsOfFace",
        (PyCFunction)create_cellsOfFace_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "create_2dfaces",
        (PyCFunction)create_2dfaces_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "create_cell_faceid",
        (PyCFunction)create_cell_faceid_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "create_3dfaces",
        (PyCFunction)create_3dfaces_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "create_NormalFacesOfCell",
        (PyCFunction)create_NormalFacesOfCell_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "face_gradient_info_2d",
        (PyCFunction)face_gradient_info_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "variables",
        (PyCFunction)variables_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "variables_3d",
        (PyCFunction)variables_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static PyModuleDef_Slot pyccel_ddm_slots[] = {
    {Py_mod_exec, exec_func},
    {0, NULL},
};

/*........................................*/

static struct PyModuleDef pyccel_ddm_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "pyccel_ddm",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    pyccel_ddm_methods,
    pyccel_ddm_slots
};

/*........................................*/

PyMODINIT_FUNC PyInit_pyccel_ddm(void)
{
    import_array();
    return PyModuleDef_Init(&pyccel_ddm_module);
}
