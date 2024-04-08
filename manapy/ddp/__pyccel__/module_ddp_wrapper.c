#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#include <stdlib.h>
#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"

void bind_c_create_info_2dfaces(int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_namen, int64_t *namen, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t nbfaces, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_mesuref, double *mesuref, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_namef, int64_t *namef);
void bind_c_create_info_3dfaces(int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_namen, int64_t *namen, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t nbfaces, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_mesuref, double *mesuref, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_namef, int64_t *namef);
void bind_c_compute_2dcentervolumeofcell(int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t nbelements, int64_t n0_center, int64_t n1_center, double *center, int64_t n0_volume, double *volume);
void bind_c_compute_3dcentervolumeofcell(int64_t n0_nodeid, int64_t n1_nodeid, int64_t *nodeid, int64_t n0_vertex, int64_t n1_vertex, double *vertex, int64_t nbelements, int64_t n0_center, int64_t n1_center, double *center, int64_t n0_volume, double *volume);
void bind_c_create_cellsofface(int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t nbelements, int64_t nbfaces, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t dim);
void bind_c_create_2dfaces(int64_t n0_nodeidc, int64_t n1_nodeidc, int64_t *nodeidc, int64_t nbelements, int64_t n0_faces, int64_t n1_faces, int64_t *faces, int64_t n0_cellf, int64_t n1_cellf, int64_t *cellf);
void bind_c_create_cell_faceid(int64_t nbelements, int64_t n0_oldTonewIndex, int64_t *oldTonewIndex, int64_t n0_cellf, int64_t n1_cellf, int64_t *cellf, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t dim);
void bind_c_create_3dfaces(int64_t n0_nodeidc, int64_t n1_nodeidc, int64_t *nodeidc, int64_t nbelements, int64_t n0_faces, int64_t n1_faces, int64_t *faces, int64_t n0_cellf, int64_t n1_cellf, int64_t *cellf);
void bind_c_create_normalfacesofcell(int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerf, int64_t n1_centerf, double *centerf, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t n0_normal, int64_t n1_normal, double *normal, int64_t nbelements, int64_t n0_nf, int64_t n1_nf, int64_t n2_nf, double *nf, int64_t dim);

/*........................................*/


/*........................................*/

/*........................................*/
PyObject *create_info_2dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t *cellid;
    int64_t *nodeid;
    int64_t *namen;
    double *vertex;
    double *centerc;
    int64_t nbfaces;
    double *normalf;
    double *mesuref;
    double *centerf;
    int64_t *namef;
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
    int64_t cellid_dim;
    int64_t cellid_dim_0001;
    int64_t nodeid_dim;
    int64_t nodeid_dim_0001;
    int64_t namen_dim;
    int64_t vertex_dim;
    int64_t vertex_dim_0001;
    int64_t centerc_dim;
    int64_t centerc_dim_0001;
    int64_t normalf_dim;
    int64_t normalf_dim_0001;
    int64_t mesuref_dim;
    int64_t centerf_dim;
    int64_t centerf_dim_0001;
    int64_t namef_dim;
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
    if (!pyarray_checker(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = PyArray_DATA(cellid_tmp);
    }
    if (!pyarray_checker(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = PyArray_DATA(nodeid_tmp);
    }
    if (!pyarray_checker(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = PyArray_DATA(namen_tmp);
    }
    if (!pyarray_checker(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = PyArray_DATA(vertex_tmp);
    }
    if (!pyarray_checker(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = PyArray_DATA(centerc_tmp);
    }
    if (PyIs_Int64(nbfaces_tmp))
    {
        nbfaces = PyInt64_to_Int64(nbfaces_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = PyArray_DATA(normalf_tmp);
    }
    if (!pyarray_checker(mesuref_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesuref = PyArray_DATA(mesuref_tmp);
    }
    if (!pyarray_checker(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = PyArray_DATA(centerf_tmp);
    }
    if (!pyarray_checker(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = PyArray_DATA(namef_tmp);
    }
    cellid_dim = PyArray_DIM(cellid_tmp, 0);
    cellid_dim_0001 = PyArray_DIM(cellid_tmp, 1);
    nodeid_dim = PyArray_DIM(nodeid_tmp, 0);
    nodeid_dim_0001 = PyArray_DIM(nodeid_tmp, 1);
    namen_dim = PyArray_DIM(namen_tmp, 0);
    vertex_dim = PyArray_DIM(vertex_tmp, 0);
    vertex_dim_0001 = PyArray_DIM(vertex_tmp, 1);
    centerc_dim = PyArray_DIM(centerc_tmp, 0);
    centerc_dim_0001 = PyArray_DIM(centerc_tmp, 1);
    normalf_dim = PyArray_DIM(normalf_tmp, 0);
    normalf_dim_0001 = PyArray_DIM(normalf_tmp, 1);
    mesuref_dim = PyArray_DIM(mesuref_tmp, 0);
    centerf_dim = PyArray_DIM(centerf_tmp, 0);
    centerf_dim_0001 = PyArray_DIM(centerf_tmp, 1);
    namef_dim = PyArray_DIM(namef_tmp, 0);
    bind_c_create_info_2dfaces(cellid_dim, cellid_dim_0001, cellid, nodeid_dim, nodeid_dim_0001, nodeid, namen_dim, namen, vertex_dim, vertex_dim_0001, vertex, centerc_dim, centerc_dim_0001, centerc, nbfaces, normalf_dim, normalf_dim_0001, normalf, mesuref_dim, mesuref, centerf_dim, centerf_dim_0001, centerf, namef_dim, namef);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_info_3dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t *cellid;
    int64_t *nodeid;
    int64_t *namen;
    double *vertex;
    double *centerc;
    int64_t nbfaces;
    double *normalf;
    double *mesuref;
    double *centerf;
    int64_t *namef;
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
    int64_t cellid_dim;
    int64_t cellid_dim_0001;
    int64_t nodeid_dim;
    int64_t nodeid_dim_0001;
    int64_t namen_dim;
    int64_t vertex_dim;
    int64_t vertex_dim_0001;
    int64_t centerc_dim;
    int64_t centerc_dim_0001;
    int64_t normalf_dim;
    int64_t normalf_dim_0001;
    int64_t mesuref_dim;
    int64_t centerf_dim;
    int64_t centerf_dim_0001;
    int64_t namef_dim;
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
    if (!pyarray_checker(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = PyArray_DATA(cellid_tmp);
    }
    if (!pyarray_checker(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = PyArray_DATA(nodeid_tmp);
    }
    if (!pyarray_checker(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = PyArray_DATA(namen_tmp);
    }
    if (!pyarray_checker(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = PyArray_DATA(vertex_tmp);
    }
    if (!pyarray_checker(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = PyArray_DATA(centerc_tmp);
    }
    if (PyIs_Int64(nbfaces_tmp))
    {
        nbfaces = PyInt64_to_Int64(nbfaces_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = PyArray_DATA(normalf_tmp);
    }
    if (!pyarray_checker(mesuref_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesuref = PyArray_DATA(mesuref_tmp);
    }
    if (!pyarray_checker(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = PyArray_DATA(centerf_tmp);
    }
    if (!pyarray_checker(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = PyArray_DATA(namef_tmp);
    }
    cellid_dim = PyArray_DIM(cellid_tmp, 0);
    cellid_dim_0001 = PyArray_DIM(cellid_tmp, 1);
    nodeid_dim = PyArray_DIM(nodeid_tmp, 0);
    nodeid_dim_0001 = PyArray_DIM(nodeid_tmp, 1);
    namen_dim = PyArray_DIM(namen_tmp, 0);
    vertex_dim = PyArray_DIM(vertex_tmp, 0);
    vertex_dim_0001 = PyArray_DIM(vertex_tmp, 1);
    centerc_dim = PyArray_DIM(centerc_tmp, 0);
    centerc_dim_0001 = PyArray_DIM(centerc_tmp, 1);
    normalf_dim = PyArray_DIM(normalf_tmp, 0);
    normalf_dim_0001 = PyArray_DIM(normalf_tmp, 1);
    mesuref_dim = PyArray_DIM(mesuref_tmp, 0);
    centerf_dim = PyArray_DIM(centerf_tmp, 0);
    centerf_dim_0001 = PyArray_DIM(centerf_tmp, 1);
    namef_dim = PyArray_DIM(namef_tmp, 0);
    bind_c_create_info_3dfaces(cellid_dim, cellid_dim_0001, cellid, nodeid_dim, nodeid_dim_0001, nodeid, namen_dim, namen, vertex_dim, vertex_dim_0001, vertex, centerc_dim, centerc_dim_0001, centerc, nbfaces, normalf_dim, normalf_dim_0001, normalf, mesuref_dim, mesuref, centerf_dim, centerf_dim_0001, centerf, namef_dim, namef);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *Compute_2dcentervolumeOfCell_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t *nodeid;
    double *vertex;
    int64_t nbelements;
    double *center;
    double *volume;
    PyArrayObject *nodeid_tmp;
    PyArrayObject *vertex_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *center_tmp;
    PyArrayObject *volume_tmp;
    int64_t nodeid_dim;
    int64_t nodeid_dim_0001;
    int64_t vertex_dim;
    int64_t vertex_dim_0001;
    int64_t center_dim;
    int64_t center_dim_0001;
    int64_t volume_dim;
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
    if (!pyarray_checker(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = PyArray_DATA(nodeid_tmp);
    }
    if (!pyarray_checker(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = PyArray_DATA(vertex_tmp);
    }
    if (PyIs_Int64(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(center_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        center = PyArray_DATA(center_tmp);
    }
    if (!pyarray_checker(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = PyArray_DATA(volume_tmp);
    }
    nodeid_dim = PyArray_DIM(nodeid_tmp, 0);
    nodeid_dim_0001 = PyArray_DIM(nodeid_tmp, 1);
    vertex_dim = PyArray_DIM(vertex_tmp, 0);
    vertex_dim_0001 = PyArray_DIM(vertex_tmp, 1);
    center_dim = PyArray_DIM(center_tmp, 0);
    center_dim_0001 = PyArray_DIM(center_tmp, 1);
    volume_dim = PyArray_DIM(volume_tmp, 0);
    bind_c_compute_2dcentervolumeofcell(nodeid_dim, nodeid_dim_0001, nodeid, vertex_dim, vertex_dim_0001, vertex, nbelements, center_dim, center_dim_0001, center, volume_dim, volume);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *Compute_3dcentervolumeOfCell_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t *nodeid;
    double *vertex;
    int64_t nbelements;
    double *center;
    double *volume;
    PyArrayObject *nodeid_tmp;
    PyArrayObject *vertex_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *center_tmp;
    PyArrayObject *volume_tmp;
    int64_t nodeid_dim;
    int64_t nodeid_dim_0001;
    int64_t vertex_dim;
    int64_t vertex_dim_0001;
    int64_t center_dim;
    int64_t center_dim_0001;
    int64_t volume_dim;
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
    if (!pyarray_checker(nodeid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeid = PyArray_DATA(nodeid_tmp);
    }
    if (!pyarray_checker(vertex_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertex = PyArray_DATA(vertex_tmp);
    }
    if (PyIs_Int64(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(center_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        center = PyArray_DATA(center_tmp);
    }
    if (!pyarray_checker(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = PyArray_DATA(volume_tmp);
    }
    nodeid_dim = PyArray_DIM(nodeid_tmp, 0);
    nodeid_dim_0001 = PyArray_DIM(nodeid_tmp, 1);
    vertex_dim = PyArray_DIM(vertex_tmp, 0);
    vertex_dim_0001 = PyArray_DIM(vertex_tmp, 1);
    center_dim = PyArray_DIM(center_tmp, 0);
    center_dim_0001 = PyArray_DIM(center_tmp, 1);
    volume_dim = PyArray_DIM(volume_tmp, 0);
    bind_c_compute_3dcentervolumeofcell(nodeid_dim, nodeid_dim_0001, nodeid, vertex_dim, vertex_dim_0001, vertex, nbelements, center_dim, center_dim_0001, center, volume_dim, volume);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_cellsOfFace_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t *faceid;
    int64_t nbelements;
    int64_t nbfaces;
    int64_t *cellid;
    int64_t dim;
    PyArrayObject *faceid_tmp;
    PyObject *nbelements_tmp;
    PyObject *nbfaces_tmp;
    PyArrayObject *cellid_tmp;
    PyObject *dim_tmp;
    int64_t faceid_dim;
    int64_t faceid_dim_0001;
    int64_t cellid_dim;
    int64_t cellid_dim_0001;
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
    if (!pyarray_checker(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = PyArray_DATA(faceid_tmp);
    }
    if (PyIs_Int64(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (PyIs_Int64(nbfaces_tmp))
    {
        nbfaces = PyInt64_to_Int64(nbfaces_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = PyArray_DATA(cellid_tmp);
    }
    if (PyIs_Int64(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    faceid_dim = PyArray_DIM(faceid_tmp, 0);
    faceid_dim_0001 = PyArray_DIM(faceid_tmp, 1);
    cellid_dim = PyArray_DIM(cellid_tmp, 0);
    cellid_dim_0001 = PyArray_DIM(cellid_tmp, 1);
    bind_c_create_cellsofface(faceid_dim, faceid_dim_0001, faceid, nbelements, nbfaces, cellid_dim, cellid_dim_0001, cellid, dim);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_2dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t *nodeidc;
    int64_t nbelements;
    int64_t *faces;
    int64_t *cellf;
    PyArrayObject *nodeidc_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *faces_tmp;
    PyArrayObject *cellf_tmp;
    int64_t nodeidc_dim;
    int64_t nodeidc_dim_0001;
    int64_t faces_dim;
    int64_t faces_dim_0001;
    int64_t cellf_dim;
    int64_t cellf_dim_0001;
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
    if (!pyarray_checker(nodeidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidc = PyArray_DATA(nodeidc_tmp);
    }
    if (PyIs_Int64(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(faces_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faces = PyArray_DATA(faces_tmp);
    }
    if (!pyarray_checker(cellf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellf = PyArray_DATA(cellf_tmp);
    }
    nodeidc_dim = PyArray_DIM(nodeidc_tmp, 0);
    nodeidc_dim_0001 = PyArray_DIM(nodeidc_tmp, 1);
    faces_dim = PyArray_DIM(faces_tmp, 0);
    faces_dim_0001 = PyArray_DIM(faces_tmp, 1);
    cellf_dim = PyArray_DIM(cellf_tmp, 0);
    cellf_dim_0001 = PyArray_DIM(cellf_tmp, 1);
    bind_c_create_2dfaces(nodeidc_dim, nodeidc_dim_0001, nodeidc, nbelements, faces_dim, faces_dim_0001, faces, cellf_dim, cellf_dim_0001, cellf);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_cell_faceid_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t nbelements;
    int64_t *oldTonewIndex;
    int64_t *cellf;
    int64_t *faceid;
    int64_t dim;
    PyObject *nbelements_tmp;
    PyArrayObject *oldTonewIndex_tmp;
    PyArrayObject *cellf_tmp;
    PyArrayObject *faceid_tmp;
    PyObject *dim_tmp;
    int64_t oldTonewIndex_dim;
    int64_t cellf_dim;
    int64_t cellf_dim_0001;
    int64_t faceid_dim;
    int64_t faceid_dim_0001;
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
    if (PyIs_Int64(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(oldTonewIndex_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldTonewIndex = PyArray_DATA(oldTonewIndex_tmp);
    }
    if (!pyarray_checker(cellf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellf = PyArray_DATA(cellf_tmp);
    }
    if (!pyarray_checker(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = PyArray_DATA(faceid_tmp);
    }
    if (PyIs_Int64(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    oldTonewIndex_dim = PyArray_DIM(oldTonewIndex_tmp, 0);
    cellf_dim = PyArray_DIM(cellf_tmp, 0);
    cellf_dim_0001 = PyArray_DIM(cellf_tmp, 1);
    faceid_dim = PyArray_DIM(faceid_tmp, 0);
    faceid_dim_0001 = PyArray_DIM(faceid_tmp, 1);
    bind_c_create_cell_faceid(nbelements, oldTonewIndex_dim, oldTonewIndex, cellf_dim, cellf_dim_0001, cellf, faceid_dim, faceid_dim_0001, faceid, dim);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_3dfaces_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    int64_t *nodeidc;
    int64_t nbelements;
    int64_t *faces;
    int64_t *cellf;
    PyArrayObject *nodeidc_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *faces_tmp;
    PyArrayObject *cellf_tmp;
    int64_t nodeidc_dim;
    int64_t nodeidc_dim_0001;
    int64_t faces_dim;
    int64_t faces_dim_0001;
    int64_t cellf_dim;
    int64_t cellf_dim_0001;
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
    if (!pyarray_checker(nodeidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidc = PyArray_DATA(nodeidc_tmp);
    }
    if (PyIs_Int64(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(faces_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faces = PyArray_DATA(faces_tmp);
    }
    if (!pyarray_checker(cellf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellf = PyArray_DATA(cellf_tmp);
    }
    nodeidc_dim = PyArray_DIM(nodeidc_tmp, 0);
    nodeidc_dim_0001 = PyArray_DIM(nodeidc_tmp, 1);
    faces_dim = PyArray_DIM(faces_tmp, 0);
    faces_dim_0001 = PyArray_DIM(faces_tmp, 1);
    cellf_dim = PyArray_DIM(cellf_tmp, 0);
    cellf_dim_0001 = PyArray_DIM(cellf_tmp, 1);
    bind_c_create_3dfaces(nodeidc_dim, nodeidc_dim_0001, nodeidc, nbelements, faces_dim, faces_dim_0001, faces, cellf_dim, cellf_dim_0001, cellf);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

/*........................................*/
PyObject *create_NormalFacesOfCell_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    double *centerc;
    double *centerf;
    int64_t *faceid;
    double *normal;
    int64_t nbelements;
    double *nf;
    int64_t dim;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerf_tmp;
    PyArrayObject *faceid_tmp;
    PyArrayObject *normal_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *nf_tmp;
    PyObject *dim_tmp;
    int64_t centerc_dim;
    int64_t centerc_dim_0001;
    int64_t centerf_dim;
    int64_t centerf_dim_0001;
    int64_t faceid_dim;
    int64_t faceid_dim_0001;
    int64_t normal_dim;
    int64_t normal_dim_0001;
    int64_t nf_dim;
    int64_t nf_dim_0001;
    int64_t nf_dim_0002;
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
    if (!pyarray_checker(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = PyArray_DATA(centerc_tmp);
    }
    if (!pyarray_checker(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = PyArray_DATA(centerf_tmp);
    }
    if (!pyarray_checker(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = PyArray_DATA(faceid_tmp);
    }
    if (!pyarray_checker(normal_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normal = PyArray_DATA(normal_tmp);
    }
    if (PyIs_Int64(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    if (!pyarray_checker(nf_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nf = PyArray_DATA(nf_tmp);
    }
    if (PyIs_Int64(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be 64 bit int\"");
        return NULL;
    }
    centerc_dim = PyArray_DIM(centerc_tmp, 0);
    centerc_dim_0001 = PyArray_DIM(centerc_tmp, 1);
    centerf_dim = PyArray_DIM(centerf_tmp, 0);
    centerf_dim_0001 = PyArray_DIM(centerf_tmp, 1);
    faceid_dim = PyArray_DIM(faceid_tmp, 0);
    faceid_dim_0001 = PyArray_DIM(faceid_tmp, 1);
    normal_dim = PyArray_DIM(normal_tmp, 0);
    normal_dim_0001 = PyArray_DIM(normal_tmp, 1);
    nf_dim = PyArray_DIM(nf_tmp, 0);
    nf_dim_0001 = PyArray_DIM(nf_tmp, 1);
    nf_dim_0002 = PyArray_DIM(nf_tmp, 2);
    bind_c_create_normalfacesofcell(centerc_dim, centerc_dim_0001, centerc, centerf_dim, centerf_dim_0001, centerf, faceid_dim, faceid_dim_0001, faceid, normal_dim, normal_dim_0001, normal, nbelements, nf_dim, nf_dim_0001, nf_dim_0002, nf, dim);
    result = Py_BuildValue("");
    return result;
}
/*........................................*/

static PyMethodDef module_ddp_methods[] = {
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
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static struct PyModuleDef module_ddp_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "module_ddp",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    -1,
    module_ddp_methods
};

/*........................................*/

PyMODINIT_FUNC PyInit_module_ddp(void)
{
    PyObject *m;
    import_array();
    m = PyModule_Create(&module_ddp_module);
    if (m == NULL) return NULL;
    return m;
}
