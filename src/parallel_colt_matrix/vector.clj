(ns parallel-colt-matrix.vector
  (:use [clojure.core.matrix.protocols])
  (:import [cern.colt.matrix AbstractMatrix1D]
           [cern.colt.matrix.tdouble.impl DenseDoubleMatrix1D]))


(extend-type AbstractMatrix1D
  PImplementation
  (implementation-key [m]
    :parallel-colt)
  (construct-matrix [m data]
    (DenseDoubleMatrix1D. (into-array Double/TYPE data)))
  (new-vector [m length]
    (DenseDoubleMatrix1D. length))
  (new-matrix [m rows columns]
    (throw (Exception. "Only vector (1D matrix) use new-vector")))
  (new-matrix-nd [m shape]
    (throw (Exception. "Only vector (1D matrix) use new-vector")))
  (supports-dimensionality? [m dimension]
    (== 1 dimension))

  PDimensionInfo
  (dimensionality [m]
    1)
  (get-shape [m]
    (.size m))
  (is-scalar? [m]
    false)
  (is-vector? [m]
    true)
  (dimension-count [m dimension-number]
    (assert (supports-dimensionality? m dimension-number))
    (get-shape m))

  PIndexedAccess
  (get-1d [m row]
    (.get m row))
  (get-2d [m row columns]
    (throw (Exception. "It is a vector with just one dimension, use get-1d")))
  (get-nd [m indexes]
    (throw (Exception. "It is a vector with just one dimension, use get-1d")))

  PIndexedSetting
  (set-1d [m row v]
    (let [other (.copy m)]
      (.set other row v)))
  (set-2d [m row col v]
    (throw (Exception. "It is a vector with just one dimension, use set-1d")))
  (set-nd [m indexes v]
    (throw (Exception. "It is a vector with just one dimension, use set-1d")))
  (is-mutable? [m]
    true)

  PIndexedSettingMutable
  (set-1d! [m row v]
    (.set m row v))
  (set-2d! [m row col v]
    (throw (Exception. "It is a vector with just one dimension, use set-1d")))
  (set-nd! [m indexes v]
    (throw (Exception. "It is a vector with just one dimension, use set-1d")))

  PMatrixCloning
  (clone [m]
    (.copy m))

  PTypeInfo
  (element-type [m]
    (class (get-1d m 1)))

  PZeroDimensionAccess
  (get-0d [m]
    (assert (= 1 (get-shape m)))
    (get-1d m 0))
  (set-0d! [m value]
    (assert (= 1 (get-shape m)))
    (set-1d! m 0 value)))

(defn get-vector
  ([] (DenseDoubleMatrix1D. 4))
  ([data] (construct-matrix (get-vector) data)))