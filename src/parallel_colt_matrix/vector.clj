(ns parallel-colt-matrix.vector
  (:use [clojure.core.matrix.protocols])
  (:import [cern.colt.matrix AbstractMatrix1D]
           [cern.colt.matrix.tdouble.impl DenseDoubleMatrix1D]
           [cern.jet.math.tdouble DoubleFunctions]
           [cern.colt.matrix.tdouble.algo DenseDoubleAlgebra]))

(defn get-vector
  ([] (DenseDoubleMatrix1D. 4))
  ([data] (construct-matrix (get-vector) data)))

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
    [(.size m)])
  (is-scalar? [m]
    false)
  (is-vector? [m]
    true)
  (dimension-count [m dimension-number]
    (assert (supports-dimensionality? m dimension-number))
    (.size m))

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
      (.set other row v)
      other))
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
    (assert (== 1 (dimension-count m 1)))
    (assert (= [1] (get-shape m)))
    (get-1d m 0))
  (set-0d! [m value]
    (assert (== 1 (dimension-count m 1)))
    (assert (= [1] (get-shape m)))
    (set-1d! m 0 value))

  PCoercion
  (coerce-param [m param]
    (let [param (flatten param)]
      (construct-matrix m param)))

  PBroadcast
  (broadcast [m target-shape]
    nil)

  PConversion
  (convert-to-nested-vectors [m]
    (vec (.elements m)))

  PSubVector
  (subvector [m start len]
    (let [v (convert-to-nested-vectors m)]
      (->> (subvec v start (+ start len))
           (construct-matrix m ))))

  PAssignment
  (assign! [m source]
    ;; auto-check in the java class
    (.assign m (into-array Double/TYPE (flatten source))))
  (assign-array!
    ([m arr]
       (assign! m  arr))
    ([m arr start length]
       (assign! m (subvec (vec (flatten arr)) start (+ start length)))))

  PDoubleArrayOutput
  (to-double-array [m]
    "Returns a double array containing the values of m in row-major order. May or may not be
     the internal double array used by m, depending on the implementation."
    (.elements (.copy m)))
  (as-double-array [m]
    "Returns the internal double array used by m. If no such array is used, returns nil.
     Provides an opportunity to avoid copying the internal array."
    (.elements m))

  PMatrixEquality
  (matrix-equals [a b]
    (cond (instance? AbstractMatrix1D)
          (.equals a b)
          (and (satisfies? PDimensionInfo b) (satisfies? PFunctionalOperations b))
          (and (= (get-shape a) (get-shape b))
               (= (element-seq a) (element-seq b)))
          :else false))

  PVectorTransform
  (vector-transform [m v])
  (vector-transform! [m v])

  PMatrixMultiply
  (element-multiply [m a]
    (let [multiplier (. DoubleFunctions mult a)
          other (.copy m)]
      (.assign other multiplier)))

  PMatrixScaling
  (scale [m a]
    (element-multiply m a))

  PMatrixMutableScaling
  (scale! [m a]
    (.assign m (. DoubleFunctions mult a)))

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

  PMatrixAddMutable
  (matrix-add! [m a]
    (assert (= (get-shape m) (get-shape a)))
    (let [sum (. DoubleFunctions plus)]
      (.assign m a sum)))
  (matrix-sub! [m a]
    (assert (= (get-shape m) (get-shape a)))
    (let [minus (. DoubleFunctions minus)]
      (.assign m a minus)))
  
  PVectorOps
  (vector-dot [a b]
     "Dot product of two vectors. Should return a scalar value."
    (.zDotProduct a b))
  (length [a]
     "Euclidian length of a vector."
    (.norm2 (DenseDoubleAlgebra.) a))
  (length-squared [a]
     "Squared Euclidean length of a vector."
    (.mult (DenseDoubleAlgebra.) a a))
  (normalise [a]
    "Returns a new vector, normalised to length 1.0"
    (let [other (.copy a)]
      (.normalize other)
      other))

  PVectorCross
  (cross-product [a b]
    "Cross product of two vectors"
    (condp = (mapv get-shape [a b])
      [[2] [2]] (let [[a1 a2] (element-seq a)
                  [b1 b2] (element-seq b)]
              (- (* a1 b2) (* a2 b1)))
      [[3] [3]] (let [[a1 a2 a3] (element-seq a)
                  [b1 b2 b3] (element-seq b)]
              (get-vector [(- (* a2 b3) (* a3 b2))
                           (- (* a3 b1) (* a1 b3))
                           (- (* a1 b2) (* a2 b1))]))
      :else (throw (Exception. "The cross product need the vector to be of the same dimension of either 2 or 3"))))
  (cross-product! [a b] ;; Possible error or not expected behaviour in case of 2dimension vector that return a single element
    "Calculate cross product of two vectors, storing the result in the first vector"
    (.assign a (cross-product a b))) 

  PMutableVectorOps
  (normalise! [a]
    (.normalize a)
    a)
  
  PFunctionalOperations
  (element-seq [m]
    (for [i (range (dimension-count m 1))]
      (get-1d m i)))

  )

;;;; add VectorOps-test
