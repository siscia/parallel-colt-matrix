(ns parallel-colt-matrix.vector-test
  (:use [clojure.test]
        [clojure.core.matrix.protocols]
        [clojure.core.matrix.compliance-tester])
  (:require [parallel-colt-matrix.vector :as pcv])
  (:import [cern.colt.matrix.tdouble.impl DenseDoubleMatrix1D]))

(def v (pcv/get-vector [1 2 3 4]))

(deftest Implementation-test
  (is (= :parallel-colt (implementation-key v)))
  (is (= v (construct-matrix v [1 2 3 4])))
  (is (supports-dimensionality? v 1))
  (is (not (supports-dimensionality? v 2))))

(deftest DimensionInfo-test
  (is (= 1 (dimensionality v)))
  (is (= 4 (get-shape v)))
  (is (not (is-scalar? v)))
  (is (is-vector? v))
  (is (= 4 (dimension-count v 1))))