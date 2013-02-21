(ns parallel-colt-matrix.core-test
  (:use clojure.test)
  (:require [parallel-colt-matrix.core :as pc])
  (:import [cern.colt.matrix.tdouble.impl DenseDoubleMatrix2D]))

(def mpc (DenseDoubleMatrix2D. 2 2))
             
(deftest vector-dimensionality-test
  (is (== 2 (pc/vector-dimensionality [[1 2 3]])))
  (is (== 3 (pc/vector-dimensionality [[[1] [2] [3] [4]]]))))

(deftest implementation-key-test
  (is (= :parallel-colt (pc/implementation-key mpc))))

;; (deftest construct-test
;;   (is (= (construct-matrix mpc [[1 2 3] [3 4 5]]) [[1 2 3] [3 4 5]])))

;; (deftest test-2-2-sparse
;;   (let [m1 (DenseDoubleMatrix2D. 2 2)]
;;     (println m1 (class m1))
;;     (compliance-test m1)))