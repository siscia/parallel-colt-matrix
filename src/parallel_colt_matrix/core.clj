(ns parallel-colt-matrix.core
  (:use [core.matrix.protocols])
  (:require [core.matrix.implementations :as imp])
  (:import [cern.colt.matrix.tdouble DoubleMatrix2D])
  (:import [cern.colt.matrix.tdouble.impl DenseDoubleMatrix2D DiagonalDoubleMatrix2D DenseDoubleMatrix1D])
  (:import [cern.colt.function.tdouble DoubleFunction DoubleDoubleFunction])
  (:import [cern.jet.math.tdouble DoubleFunctions])
  (:import [cern.colt.matrix.tdouble.algo DenseDoubleAlgebra]))

(defn vector-dimensionality ^long [m]
  "Calculates the dimensionality (== nesting depth) of nested persistent vectors"
  (cond
    (vector? m)
    (loop [m m c 0]
      (if (vector? m)
        (recur (first m) (inc c))
        c))
    (satisfies? PDimensionInfo m) (long (dimensionality m))
    :else (throw (Exception. (str "Don't know how to find dimension, " (class m) "is not a vactor nor it implement PDimensionInfo")))))

(extend DoubleMatrix2D
  PImplementation
  {:implementation-key (fn [_] :parallel-colt)})

(extend-type DoubleMatrix2D
  PImplementation
  (implementation-key [m]
    "Returns a keyword representing this implementation. 
     Each implementation should have one unique key."
    :parallel-colt)
  (construct-matrix [m data]
    "Returns a new matrix containing the given data. Data should be in the form of either nested sequences or a valid existing matrix"
    (cond
     (isa? DoubleMatrix2D data) data
     (and (vector? data) (== 2 (vector-dimensionality data))) (DenseDoubleMatrix2D. (into-array (map #(into-array Double/TYPE (map double %)) data)))
     (satisfies? PConversion data) (construct-matrix m (convert-to-nested-vectors data))
     :else (throw (Exception. (str "Don't know how to convert " (class data) "into a 2D vector, it need to implement the PCoversion.")))))
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
   (supports-dimensionality? [m dimensions]
     "Returns true if the implementation supports matrices with the given number of dimensions."
     (== 2 dimensions))

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

  PIndexedSetting
  (set-1d [m row v]
    (throw (Exception. "It is a 2D matrix, specify another dimension, use set-2d(!)")))
  (set-2d [m row column v]
    (let [other (.copy m)]  ;;Keeping immutability
      (.set other row column v)
      other)) ;;We could use .setQuick but then we need
  ;;to add a pre-condition to check if the index is in the bound, and
  ;;I guess it will be slower than simply use .set
  (set-nd [m indexes v]
    (throw (Exception. "It is a 2D matrix, no other dimension, use set-2d(!)")))
  (is-mutable? [m]
    true)
  
  PSpecialisedConstructors
  (identity-matrix [m dims] "Create a 2D identity matrix with the given number of dimensions"
    (DenseDoubleMatrix2D. dims dims)) ;;Square matrix
  (diagonal-matrix [m diagonal-values] "Create a diagonal matrix with the specified leading diagonal values"
    (let [dim (count diagonal-values)
          m (DiagonalDoubleMatrix2D. dim dim 0)]
      (dotimes [i dim]
        (.setQuick m i i (double (nth diagonal-values i))))
      m))

  PCoercion
  (coerce-param [m param]
    "Attempts to coerce param into a matrix format supported by the implementation of matrix m.
     May return nil if unable to do so, in which case a default implementation can be used."
    (condp == (dimensionality param)
      0 param
      1 (DenseDoubleMatrix1D. (double-array param))  ;;TODO wait to
      ;;implement also 1D Matrix in PColt, good enough for now
      2 (construct-matrix m param)
      (throw (Exception. "Need to be done"))))
  
  PMatrixEquality
  (matrix-equals [a b]
    (.equals a b))
  
  ;; PAssignment ;;TODO
  ;; "Protocol for assigning values to mutable matrices."
  ;; (assign-array! [m arr] "Sets all the values in a matrix from a Java array, in row-major order")
  ;; (assign! [m source] "Sets all the values in a matrix from a matrix source")

  PMatrixMultiply
  (matrix-multiply [m a]
    (.zMult m (coerce-param m a) nil))
  (element-multiply [m a]
    (let [multiplier (. DoubleFunctions mult a)
          other (.copy m)]
      (.assign other multiplier)))

  ;; PVectorTransform ;;TODO
  ;; "Protocol to support transformation of a vector to another vector. 
  ;;  Is equivalent to matrix multiplication when 2D matrices are used as transformations.
  ;;  But other transformations are possible, e.g. affine transformations."
  ;; (vector-transform [m v] "Transforms a vector")
  ;; (vector-transform! [m v] "Transforms a vector in place - mutates the vector argument")

  ;; PMatrixScaling ;;TODO
  ;; "Protocol to support matrix scaling by scalar values"
  ;; (scale [m a])
  ;; (pre-scale [m a])
  
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

  ;; PVectorOps ;;TODO
  ;; "Protocol to support common vector operations."
  ;; (vector-dot [a b])
  ;; (length [a])
  ;; (length-squared [a])
  ;; (normalise [a])

  PMatrixOps
  (trace [m]
    (.trace (DenseDoubleAlgebra.) m))
  (determinant [m]
    (.det (DenseDoubleAlgebra.) m))
  (inverse [m]
    (.inverse (DenseDoubleAlgebra.) m))
  (negate [m]
    (element-multiply m -1))
  (transpose [m]
    (.transpose (DenseDoubleAlgebra.) m))

  ;; PMathsFunctions ;;TODO

  PMatrixSlices ;;TODO
  (get-row [m i]
    (-> (.viewRow m i) (.toArray) (vec))) ;; same as get-1d
  (get-column [m i]
    (-> (.viewColumn m i) (.toArray) (vec)))
  (get-major-slice [m i]
    (get-row m i)) ;; ???? DON'T GET IT
  (get-slice [m dimension i]
    (condp == dimension
      0 (get-row m i)
      1 (get-column m i)
      (throw (Exception. "It is a 2D matrix, it only have 2 dimension"))))
  
  PFunctionalOperations  ;;TODO
  ;; "Protocol to allow functional-style operations on matrix elements."
  ;; ;; note that protocols don't like variadic args, so we convert to
  ;; regularargs 
  (element-seq [m]
  (let [shape (get-shape m)]
    (for [row (range (shape 0))
          col (range (shape 1))]
      (get-2d m row col))))
  (element-map
    ([m f]
       (let [fun (reify DoubleFunction
                   (apply [m n] (f n)))
             other (.copy m)]
         (.assign other fun))))
  (element-map
    ([m f a]
       (let [fun (reify DoubleDoubleFunction
                   (apply [m n a] (f n a)))
             other (.copy m)]
         (.assign other fun))))
    
  ;; (element-map [m f]
  ;;              [m f a]
  ;;              [m f a more])
  ;; (element-map! [m f]
  ;;               [m f a]
  ;;               [m f a more])
  ;; (element-reduce [m f] [m f init])

  PConversion
  (convert-to-nested-vectors [m]
    (->> (.toArray m) (map vec) vec)))

;; (defn element-s "just an idea, it does not work, it might be useful
;;   for high dimension matrix"
;;   [m]
;;   (let [shape (get-shape m)]
;;     (for [[row col] (map #(range %) (get-shape m))]
;;       (get-2d m row col))))

;; (defn map-over-nested [f ar]
;;   (cond (> (vector-dimensionality ar) 1)
;;           (map #(map-over-nested f %) ar)
;;         (== 1 (vector-dimensionality ar))
;;           (map f ar)))

(imp/register-implementation (DenseDoubleMatrix2D. 2 2))
