#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include <stdlib.h>
#include "ndarrays.h"
#include "cwrapper_ndarrays.h"
#include <stdint.h>


void bind_c_initialisation_gaussian_2d(int64_t n0_ne, double *ne, int64_t n0_u, double *u, int64_t n0_v, double *v, int64_t n0_w, double *w, int64_t n0_P, double *P, int64_t n0_center, int64_t n1_center, double *center, double Pinit);
void bind_c_initialisation_gaussian_3d(int64_t n0_ne, double *ne, int64_t n0_u, double *u, int64_t n0_v, double *v, int64_t n0_w, double *w, int64_t n0_P, double *P, int64_t n0_center, int64_t n1_center, double *center, double Pinit);
double bind_c_time_step(int64_t n0_u, double *u, int64_t n0_v, double *v, int64_t n0_w, double *w, double cfl, int64_t n0_normal, int64_t n1_normal, double *normal, int64_t n0_mesure, double *mesure, int64_t n0_volume, double *volume, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t dim, double Dxx, double Dyy, double Dzz);
void bind_c_update_new_value(int64_t n0_ne_c, double *ne_c, int64_t n0_u_c, double *u_c, int64_t n0_v_c, double *v_c, int64_t n0_P_c, double *P_c, int64_t n0_rez_ne, double *rez_ne, int64_t n0_dissip_ne, double *dissip_ne, int64_t n0_src_ne, double *src_ne, double dtime, int64_t n0_vol, double *vol);

/*........................................*/


/*........................................*/

/*........................................*/
PyObject *initialisation_gaussian_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray ne = {.shape = NULL};
    t_ndarray u = {.shape = NULL};
    t_ndarray v = {.shape = NULL};
    t_ndarray w = {.shape = NULL};
    t_ndarray P = {.shape = NULL};
    t_ndarray center = {.shape = NULL};
    double Pinit;
    PyArrayObject *ne_tmp;
    PyArrayObject *u_tmp;
    PyArrayObject *v_tmp;
    PyArrayObject *w_tmp;
    PyArrayObject *P_tmp;
    PyArrayObject *center_tmp;
    PyObject *Pinit_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "ne",
        "u",
        "v",
        "w",
        "P",
        "center",
        "Pinit",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O", kwlist, &PyArray_Type, &ne_tmp, &PyArray_Type, &u_tmp, &PyArray_Type, &v_tmp, &PyArray_Type, &w_tmp, &PyArray_Type, &P_tmp, &PyArray_Type, &center_tmp, &Pinit_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(ne_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        ne = pyarray_to_ndarray(ne_tmp);
    }
    if (!pyarray_check(u_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u = pyarray_to_ndarray(u_tmp);
    }
    if (!pyarray_check(v_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v = pyarray_to_ndarray(v_tmp);
    }
    if (!pyarray_check(w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w = pyarray_to_ndarray(w_tmp);
    }
    if (!pyarray_check(P_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        P = pyarray_to_ndarray(P_tmp);
    }
    if (!pyarray_check(center_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        center = pyarray_to_ndarray(center_tmp);
    }
    if (PyIs_NativeFloat(Pinit_tmp))
    {
        Pinit = PyDouble_to_Double(Pinit_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    bind_c_initialisation_gaussian_2d(nd_ndim(&ne, 0), nd_data(&ne), nd_ndim(&u, 0), nd_data(&u), nd_ndim(&v, 0), nd_data(&v), nd_ndim(&w, 0), nd_data(&w), nd_ndim(&P, 0), nd_data(&P), nd_ndim(&center, 0), nd_ndim(&center, 1), nd_data(&center), Pinit);
    result = Py_BuildValue("");
    free_pointer(ne);
    free_pointer(u);
    free_pointer(v);
    free_pointer(w);
    free_pointer(P);
    free_pointer(center);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *initialisation_gaussian_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray ne = {.shape = NULL};
    t_ndarray u = {.shape = NULL};
    t_ndarray v = {.shape = NULL};
    t_ndarray w = {.shape = NULL};
    t_ndarray P = {.shape = NULL};
    t_ndarray center = {.shape = NULL};
    double Pinit;
    PyArrayObject *ne_tmp;
    PyArrayObject *u_tmp;
    PyArrayObject *v_tmp;
    PyArrayObject *w_tmp;
    PyArrayObject *P_tmp;
    PyArrayObject *center_tmp;
    PyObject *Pinit_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "ne",
        "u",
        "v",
        "w",
        "P",
        "center",
        "Pinit",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O", kwlist, &PyArray_Type, &ne_tmp, &PyArray_Type, &u_tmp, &PyArray_Type, &v_tmp, &PyArray_Type, &w_tmp, &PyArray_Type, &P_tmp, &PyArray_Type, &center_tmp, &Pinit_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(ne_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        ne = pyarray_to_ndarray(ne_tmp);
    }
    if (!pyarray_check(u_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u = pyarray_to_ndarray(u_tmp);
    }
    if (!pyarray_check(v_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v = pyarray_to_ndarray(v_tmp);
    }
    if (!pyarray_check(w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w = pyarray_to_ndarray(w_tmp);
    }
    if (!pyarray_check(P_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        P = pyarray_to_ndarray(P_tmp);
    }
    if (!pyarray_check(center_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        center = pyarray_to_ndarray(center_tmp);
    }
    if (PyIs_NativeFloat(Pinit_tmp))
    {
        Pinit = PyDouble_to_Double(Pinit_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native float\"");
        return NULL;
    }
    bind_c_initialisation_gaussian_3d(nd_ndim(&ne, 0), nd_data(&ne), nd_ndim(&u, 0), nd_data(&u), nd_ndim(&v, 0), nd_data(&v), nd_ndim(&w, 0), nd_data(&w), nd_ndim(&P, 0), nd_data(&P), nd_ndim(&center, 0), nd_ndim(&center, 1), nd_data(&center), Pinit);
    result = Py_BuildValue("");
    free_pointer(ne);
    free_pointer(u);
    free_pointer(v);
    free_pointer(w);
    free_pointer(P);
    free_pointer(center);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *time_step_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray u = {.shape = NULL};
    t_ndarray v = {.shape = NULL};
    t_ndarray w = {.shape = NULL};
    double cfl;
    t_ndarray normal = {.shape = NULL};
    t_ndarray mesure = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray faceid = {.shape = NULL};
    int64_t dim;
    double Dxx;
    double Dyy;
    double Dzz;
    double dt;
    PyArrayObject *u_tmp;
    PyArrayObject *v_tmp;
    PyArrayObject *w_tmp;
    PyObject *cfl_tmp;
    PyArrayObject *normal_tmp;
    PyArrayObject *mesure_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *faceid_tmp;
    PyObject *dim_tmp;
    PyObject *Dxx_tmp;
    PyObject *Dyy_tmp;
    PyObject *Dzz_tmp;
    PyObject *dt_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "u",
        "v",
        "w",
        "cfl",
        "normal",
        "mesure",
        "volume",
        "faceid",
        "dim",
        "Dxx",
        "Dyy",
        "Dzz",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!OO!O!O!O!OOOO", kwlist, &PyArray_Type, &u_tmp, &PyArray_Type, &v_tmp, &PyArray_Type, &w_tmp, &cfl_tmp, &PyArray_Type, &normal_tmp, &PyArray_Type, &mesure_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &faceid_tmp, &dim_tmp, &Dxx_tmp, &Dyy_tmp, &Dzz_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(u_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u = pyarray_to_ndarray(u_tmp);
    }
    if (!pyarray_check(v_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v = pyarray_to_ndarray(v_tmp);
    }
    if (!pyarray_check(w_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w = pyarray_to_ndarray(w_tmp);
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
    if (PyIs_NativeInt(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
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
    dt = bind_c_time_step(nd_ndim(&u, 0), nd_data(&u), nd_ndim(&v, 0), nd_data(&v), nd_ndim(&w, 0), nd_data(&w), cfl, nd_ndim(&normal, 0), nd_ndim(&normal, 1), nd_data(&normal), nd_ndim(&mesure, 0), nd_data(&mesure), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid), dim, Dxx, Dyy, Dzz);
    dt_tmp = Double_to_NumpyDouble(&dt);
    result = Py_BuildValue("O", dt_tmp);
    Py_DECREF(dt_tmp);
    free_pointer(u);
    free_pointer(v);
    free_pointer(w);
    free_pointer(normal);
    free_pointer(mesure);
    free_pointer(volume);
    free_pointer(faceid);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *update_new_value_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray ne_c = {.shape = NULL};
    t_ndarray u_c = {.shape = NULL};
    t_ndarray v_c = {.shape = NULL};
    t_ndarray P_c = {.shape = NULL};
    t_ndarray rez_ne = {.shape = NULL};
    t_ndarray dissip_ne = {.shape = NULL};
    t_ndarray src_ne = {.shape = NULL};
    double dtime;
    t_ndarray vol = {.shape = NULL};
    PyArrayObject *ne_c_tmp;
    PyArrayObject *u_c_tmp;
    PyArrayObject *v_c_tmp;
    PyArrayObject *P_c_tmp;
    PyArrayObject *rez_ne_tmp;
    PyArrayObject *dissip_ne_tmp;
    PyArrayObject *src_ne_tmp;
    PyObject *dtime_tmp;
    PyArrayObject *vol_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "ne_c",
        "u_c",
        "v_c",
        "P_c",
        "rez_ne",
        "dissip_ne",
        "src_ne",
        "dtime",
        "vol",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!OO!", kwlist, &PyArray_Type, &ne_c_tmp, &PyArray_Type, &u_c_tmp, &PyArray_Type, &v_c_tmp, &PyArray_Type, &P_c_tmp, &PyArray_Type, &rez_ne_tmp, &PyArray_Type, &dissip_ne_tmp, &PyArray_Type, &src_ne_tmp, &dtime_tmp, &PyArray_Type, &vol_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(ne_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        ne_c = pyarray_to_ndarray(ne_c_tmp);
    }
    if (!pyarray_check(u_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u_c = pyarray_to_ndarray(u_c_tmp);
    }
    if (!pyarray_check(v_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_c = pyarray_to_ndarray(v_c_tmp);
    }
    if (!pyarray_check(P_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        P_c = pyarray_to_ndarray(P_c_tmp);
    }
    if (!pyarray_check(rez_ne_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rez_ne = pyarray_to_ndarray(rez_ne_tmp);
    }
    if (!pyarray_check(dissip_ne_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dissip_ne = pyarray_to_ndarray(dissip_ne_tmp);
    }
    if (!pyarray_check(src_ne_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        src_ne = pyarray_to_ndarray(src_ne_tmp);
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
    bind_c_update_new_value(nd_ndim(&ne_c, 0), nd_data(&ne_c), nd_ndim(&u_c, 0), nd_data(&u_c), nd_ndim(&v_c, 0), nd_data(&v_c), nd_ndim(&P_c, 0), nd_data(&P_c), nd_ndim(&rez_ne, 0), nd_data(&rez_ne), nd_ndim(&dissip_ne, 0), nd_data(&dissip_ne), nd_ndim(&src_ne, 0), nd_data(&src_ne), dtime, nd_ndim(&vol, 0), nd_data(&vol));
    result = Py_BuildValue("");
    free_pointer(ne_c);
    free_pointer(u_c);
    free_pointer(v_c);
    free_pointer(P_c);
    free_pointer(rez_ne);
    free_pointer(dissip_ne);
    free_pointer(src_ne);
    free_pointer(vol);
    return result;
}
/*........................................*/

static int exec_func(PyObject* m)
{
    return 0;
}

/*........................................*/

static PyMethodDef pyccel_tools_methods[] = {
    {
        "initialisation_gaussian_2d",
        (PyCFunction)initialisation_gaussian_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "initialisation_gaussian_3d",
        (PyCFunction)initialisation_gaussian_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "time_step",
        (PyCFunction)time_step_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "update_new_value",
        (PyCFunction)update_new_value_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static PyModuleDef_Slot pyccel_tools_slots[] = {
    {Py_mod_exec, exec_func},
    {0, NULL},
};

/*........................................*/

static struct PyModuleDef pyccel_tools_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "pyccel_tools",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    pyccel_tools_methods,
    pyccel_tools_slots
};

/*........................................*/

PyMODINIT_FUNC PyInit_pyccel_tools(void)
{
    import_array();
    return PyModuleDef_Init(&pyccel_tools_module);
}
