(ns parallel-colt-matrix.vector-test
  (:use [clojure.test]
        [clojure.core.matrix.protocols]
        [clojure.core.matrix.compliance-tester])
  (:require [parallel-colt-matrix.vector :as pcv])
  (:import [cern.colt.matrix.tdouble.impl DenseDoubleMatrix1D]))

(def values [1.0 2.0 3.0 4.0])
(def v (pcv/get-vector values))

(deftest Implementation-test
  (is (= :parallel-colt (implementation-key v)))
  (is (= v (construct-matrix v values)))
  (is (supports-dimensionality? v 1))
  (is (not (supports-dimensionality? v 2))))

(deftest DimensionInfo-test
  (is (= 1 (dimensionality v)))
  (is (= [4] (get-shape v)))
  (is (not (is-scalar? v)))
  (is (is-vector? v))
  (is (= 4 (dimension-count v 1))))

(deftest IndexAccess-test
  (is (= 1.0 (get-1d v 0)))
  (is (= 3.0 (get-1d v 2)))
  (is (thrown? Exception (get-2d v 2 3)))
  (is (thrown? Exception (get-nd v [2 3 4]))))

(deftest IndexSetting-test
  (is (= 42.0 (-> (set-1d v 1 42)
                  (get-1d 1))))
  (is (= 39.0 (get-1d (set-1d v 2 39.0) 2)))
  (is (thrown? Exception (set-2d v 2 3 2)))
  (is (thrown? Exception (set-nd v [2 3 4] 4)))
  (is (thrown? Exception (set-1d v 3 'foo)))
  (is (thrown? Exception (set-2d v 2 3 'ciao)))
  (is (thrown? Exception (set-nd v [2 3 4] 'hello))))

(deftest IndexSettingMutable-test
  (let [vpc (pcv/get-vector [1 2 3 4 5 6])]
    (set-1d! vpc 0 23)
    (set-1d! vpc 1 32)
    (is (not= 1.0 (get-1d vpc 0)))
    (is (== 23.0 (get-1d vpc 0)))
    (is (not= 2.0 (get-1d vpc 1)))
    (is (== 32.0 (get-1d vpc 1)))
    (is (thrown? Exception (set-2d! v 2 3 2)))
    (is (thrown? Exception (set-nd! v [2 3 4] 4)))
    (is (thrown? Exception (set-1d! v 3 'foo)))
    (is (thrown? Exception (set-2d! v 2 3 'ciao)))
    (is (thrown? Exception (set-nd! v [2 3 4] 'hello)))))

(deftest MatrixCloning-test
  (let [vpc (pcv/get-vector [1 2 3 4])]
    (is (= vpc (clone vpc)))
    (is (not= vpc (set-1d! (clone vpc) 2 3)))))

(deftest TypeInfo-test
  (let [ex [1 2 3]
        vpc (pcv/get-vector ex)]
    (is (= java.lang.Double (element-type vpc)))))

(deftest ZeroDimensionAccess-test
  (let [vpc (pcv/get-vector [1])]
    (is (= 1.0 (get-0d vpc)))
    (set-0d! vpc 2)
    (is (= 2.0 (get-0d vpc)))
    (is (thrown? Exception (set-0d! vpc 'ciao)))
    (is (thrown? AssertionError (get-0d v)))
    (is (thrown? AssertionError (set-0d! v 3)))))

(deftest Coercion-test
  (is (= v (coerce-param v [[1 [2] 3] [[4]]]))))

(deftest Conversion-test
  (is (= values (convert-to-nested-vectors v)))
  (is (= [2.0 6.0 4.0] (-> (pcv/get-vector [2 6 4])
                           convert-to-nested-vectors))))

(deftest SubVector-test
  (is (= (pcv/get-vector [2 3]) (subvector v 1 2)))
  (let [vpc (pcv/get-vector (range 10))]
    (is (= (pcv/get-vector [3 4 5 6]) (subvector vpc 3 4)))
    (is (= (pcv/get-vector [5 6 7]) (subvector vpc 5 3))))
  (is (thrown? Exception (subvector v 5 1)))
  (is (thrown? Exception (subvector v 2 6))))

(deftest Assigment-test
  (let [vpc (pcv/get-vector [1 2 3])]
    (is (= (pcv/get-vector [6 7 8]) (assign! vpc [6 7 8])))
    (is (= (pcv/get-vector [10 11 12]) (assign! vpc [[10 [11] 12]])))
    (is (not= vpc (pcv/get-vector [1 2 3])))
    (is (thrown? IllegalArgumentException (assign! vpc [1 2])))
    (is (thrown? IllegalArgumentException (assign! vpc [1 2 3 4])))
    (is (= (pcv/get-vector [1 2 3]) (assign-array! vpc [[1 2 3]])))
    (is (= (pcv/get-vector [3 4 5]) (assign-array! vpc [[1 2] [3 4] [5 6]] 2 3)))))

(deftest DoubleArrayOutput-test
  (let [vpc (pcv/get-vector [1 2 3])]
    (is (= (vec (to-double-array vpc)) (vec (as-double-array vpc))))
    (is (= (vec (to-double-array vpc)) (vec (double-array [1 2 3]))))
    (aset (as-double-array vpc) 1 42.0)
    (is (= (pcv/get-vector [1 42 3]) vpc))
    (aset (to-double-array vpc) 1 24.0)
    (is (= (pcv/get-vector [1 42 3]) vpc))))

(deftest MatrixEquality-test
  (is (= (pcv/get-vector [1 2 3]) (pcv/get-vector [1 2 3])))
  (is (= (construct-matrix v [2 3 4]) (pcv/get-vector [2 3 4])))
  (is (matrix-equals (pcv/get-vector [1 2 3]) [1.0 2.0 3.0]))
  (is (matrix-equals (pcv/get-vector [1 2 3]) (pcv/get-vector [1 2 3])))
  (is (matrix-equals (construct-matrix v [2 3 4]) (pcv/get-vector [2 3 4]))))

(deftest MatrixScaling-test
  (let [h (pcv/get-vector [1 2 3])]
    (is (= (scale h 2) (pcv/get-vector [2 4 6])))
    (is (= (scale h -1) (pcv/get-vector [-1 -2 -3])))))

(deftest MatrixMutableScaling-test
  (let [h (pcv/get-vector [1 2 3])]
    (scale! h 2)
    (is (= h (pcv/get-vector [2 4 6])))
    (is (not= h (pcv/get-vector [1 2 3])))
    (is (= (scale! h -2) (pcv/get-vector [-4 -8 -12])))))

(deftest VectorCross-test
  (is (= -1.0 (cross-product (pcv/get-vector [1 2]) (pcv/get-vector [2 3]))))
  (is (= -251.0 (cross-product (pcv/get-vector [1 5]) (pcv/get-vector [43 -36]))))
  (is (= (pcv/get-vector [-3 6 -3])
         (cross-product (pcv/get-vector [1 2 3])
                        (pcv/get-vector [4 5 6]))))
  (is (= (pcv/get-vector [4 -8 4])
         (cross-product (pcv/get-vector [3 2 1])
                        (pcv/get-vector [1 2 3]))))
  (is (thrown? Exception (cross-product
                          (pcv/get-vector [1 2 3 4])
                          (pcv/get-vector [1 2 3 4]))))
  (is (thrown? Exception (cross-product
                          (pcv/get-vector [1 2])
                          (pcv/get-vector [1 2 3]))))
  (is (thrown? Exception (cross-product
                          (pcv/get-vector [1 2 4])
                          (pcv/get-vector [1]))))
  (let [a (pcv/get-vector [1 2])
        b (pcv/get-vector [1 2 3])]
    (cross-product! a (pcv/get-vector [2 3]))
    (is (= a (pcv/get-vector [-1 -1]))) ;;Not sure if this is right, it has to be either documented or deleted
    (cross-product! b (pcv/get-vector [9 8 7]))
    (is (= b (pcv/get-vector [-10 20 -10])))))