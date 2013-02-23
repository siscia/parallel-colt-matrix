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

  ;; (deftest test-2-2-sparse
  ;;   (compliance-test mpc)))