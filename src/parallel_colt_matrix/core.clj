(ns parallel-colt-matrix.core
  (:use [core.matrix.protocols])
  (:require [core.matrix.implementations :as imp])
  (:import [cern.colt.matrix.tdouble.impl DenseDoubleMatrix2D])
  (:import [cern.jet.math.tdouble DoubleFunctions]))


(extend-type DenseDoubleMatrix2D
  PImplementation
  (implementation-key [m]
    "Returns a keyword representing this implementation. 
     Each implementation should have one unique key."
    :parallel-colt)
  (construct-matrix [m data]
    "Returns a new matrix containing the given data. Data should be in the form of either nested sequences or a valid existing matrix"
    (DenseDoubleMatrix2D. (into-array (map #(into-array Double/TYPE (map double %)) data))))
  (new-vector [m length]
    "Returns a new vector (1D column matrix) of the given length."
    (throw (Exception. "Only 2D matrix, use new-matrix")))
  (new-matrix [m rows columns]
    "Returns a new matrix (regular 2D matrix) with the given number of rows and columns."
    (DenseDoubleMatrix2D. rows columns))
  (new-matrix-nd [m shape]
    "Returns a new general matrix of the given shape.
     Shape must be a sequence of dimension sizes."
    (throw (Exception. "Only 2D matrix, use new-matrix")))

  PDimensionInfo
  (dimensionality [m] 
    "Returns the number of dimensions of a matrix"
    2)
  (get-shape [m]
    "Returns the shape of the matrix, as an array or sequence of dimension sizes"
    [(.rows m) (.columns m)])
  (is-scalar? [m] 
    "Tests whether an object is a scalar value"
    false)
  (is-vector? [m] 
    "Tests whether an object is a vector (1D matrix)"
    false)
  (dimension-count [m dimension-number]
    "Returns the size of a specific dimension
I assumed 0 for colunms 1 for rows"
    (assert (< dimension-number 3))
    (assert (>= dimension-number 0))
    (case dimension-number
      0 (.rows m)
      1 (.columns m)))

  PIndexedAccess
  (get-1d [m row]
    (-> (.viewRow m row) (.toArray) (vec)))
  (get-2d [m row column]
    (.get m row column))
  (get-nd [m indexes]
    (throw (Exception. "It is only a 2D Array")))

  PCoercion
  (coerce-param [m param]
    "Attempts to coerce param into a matrix format supported by the implementation of matrix m.
     May return nil if unable to do so, in which case a default implementation can be used."
    nil)

  PIndexedSetting
  (set-1d [m row v]
    (throw (Exception. "It is a 2D matrix, specify another dimension, use set-2d(!)")))
  (set-2d [m row column v]
    (.set m row column v)) ;;We could use .setQuick but then we need
  ;;to add a pre-condition to check if the index is in the bound, and
  ;;I guess it will be slower than simply use .set
  (set-nd [m indexes v]
    (throw (Exception. "It is a 2D matrix, no other dimension, use set-2d(!)")))

  PMatrixEquality
  (matrix-equals [a b]
    (.equals a b))

  PMatrixAdd
  (matrix-add [m a]
    (assert (= (get-shape m) (get-shape a)))
    (let [sum (. DoubleFunctions plus)
          other (.copy m)]
      (.assign other a sum)))

  (matrix-sub [m a]
    (assert (= (get-shape m) (get-shape a)))
    (let [minus (. DoubleFunctions minus)
          other (.copy m)]
      (.assign other a minus)))

  PMatrixMultiply
  (matrix-multiply [m a]
    (.zMult m a nil))
  (element-multiply [m a]
    (let [multiplier (. DoubleFunctions mult a)
          other (.copy m)]
      (.assign other multiplier))))

(imp/register-implementation (DenseDoubleMatrix2D. 2 2))
