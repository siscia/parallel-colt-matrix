(ns parallel-colt-matrix.core-test
  (:use [clojure.test]
        [clojure.core.matrix.protocols]
        [clojure.core.matrix.compliance-tester])
  (:require [parallel-colt-matrix.core :as pc])
  (:import [cern.colt.matrix.tdouble.impl DenseDoubleMatrix2D]))

(def mpc (DenseDoubleMatrix2D. 2 2))
(def m (construct-matrix mpc [[1 2 3] [4 5 6]]))

(deftest vector-dimensionality-test
  (is (== 2 (pc/vector-dimensionality [[1 2 3]])))
  (is (== 3 (pc/vector-dimensionality [[[1] [2] [3] [4]]]))))

(deftest implementation-key-test
  (is (= :parallel-colt (implementation-key mpc))))

(deftest construct-test
  (is (matrix-equals (construct-matrix mpc [[1 2 3] [3 4 5]]) (construct-matrix mpc  [[1 2 3] [3 4 5]]))))

(deftest sequence-test
  (is (= '(1.0 2.0 3.0 4.0 5.0 6.0) (element-seq m)))
  (is (= '(1.0 2.0 3.0) (element-seq (construct-matrix mpc [[1] [2] [3]]))))
  (is (= '(1.0 2.0 3.0) (element-seq (construct-matrix mpc [[1 2 3]])))))

(deftest map-test
  (is (= (construct-matrix mpc [[2] [3] [4]]) (element-map (construct-matrix mpc [[1] [2] [3]]) inc)))
  (is (= (construct-matrix mpc [[11] [22] [33]]) (element-map (construct-matrix mpc [[1] [2] [3]]) #(* 11 %))))
  (is (= (construct-matrix mpc [[11] [22] [33]]) (element-map (construct-matrix mpc [[1] [2] [3]]) + (construct-matrix mpc [[10] [20] [30]])))))

(deftest ZeroDimensionAccess-test
  (let [m (construct-matrix mpc [[1]])]
    (is (== (get-0d m) 1.0))
    (is (== (get-0d (set-0d! m 3)) 3.0))
    (is (== (get-0d m) 3.0))))

(deftest MatrixSubComponents
  (is (= (main-diagonal (pc/get-matrix [[1 2]
                                        [3 4]])) '(1.0 4.0)))
  (is (= (main-diagonal (pc/get-matrix [[1 2 3]
                                        [4 5 6]])) '(1.0 5.0)))
  (is (= (main-diagonal (pc/get-matrix [[1 2]
                                        [3 4]
                                        [5 6]])) '(1.0 4.0)))
  (is (= (main-diagonal (pc/get-matrix [[1 2 3]
                                        [4 5 6]
                                        [7 8 9]])) '(1.0 5.0 9.0))))

(deftest MatrixScaling-test
  (is (= (scale (pc/get-matrix [[1 2] [3 4]]) 3) (pc/get-matrix [[3 6] [9 12]])))
  (is (= (scale (pc/get-matrix [[1 2] [3 4]]) -1) (pc/get-matrix [[-1 -2] [-3 -4]])))
  (is (= (scale (pc/get-matrix [[1 2] [3 4]]) 1.5) (pc/get-matrix [[1.5 3] [4.5 6]]))))

(deftest MatrixMutableScaling-test
  (let [m (pc/get-matrix [[1 2] [3 4]])]
    (is (= (pc/get-matrix [[2 4] [6 8]]) (scale! m 2)))
    (is (= (pc/get-matrix [[2 4] [6 8]]) m))
    (is (= (pc/get-matrix [[1 2] [3 4]]) (scale! m 1/2)))
    (is (= (pc/get-matrix [[1 2] [3 4]]) m))
    (is (= (scale m 3) (scale! m 3)) "Testing consistency between scale and scale!")
    (is (not= (scale (pc/get-matrix [[1 2] [3 4]]) 2) (scale m 2)) "m is now different from the start")
    (is (not= (pc/get-matrix [[1 2] [3 4]]) m))))

(deftest Reshaping-test
  (let [m (pc/get-matrix [[1 2 3] [4 5 6]])]
    (is (= (reshape m [2 2]) (pc/get-matrix [[1 2]
                                             [3 4]])))
    (is (= (reshape m [2 3]) (pc/get-matrix [[1 2]
                                             [3 4]
                                             [5 6]])))
    (is (thrown? AssertionError (reshape m [2 4])))))

(deftest DoubleArrayOutput-test
  (let [m (pc/get-matrix [[1 2 3] [4 5 6]])
        ar (double-array (flatten [[1 2 3] [4 5 6]]))]
    (is (= (seq ar) (seq (as-double-array m))))
    (is (= (seq ar) (seq (to-double-array m))))
    (aset (as-double-array m) 1 42.0)
    (is (= m (pc/get-matrix [[1 42 3] [4 5 6]])) "Check if it really change the underneath array")
    (aset (to-double-array m) 2 53.0)
    (is (not= m (pc/get-matrix [[1 42 53] [4 5 6]])) "Check correct use of TO-double-array")
    (is (= m (pc/get-matrix [[1 42 3] [4 5 6]])))))