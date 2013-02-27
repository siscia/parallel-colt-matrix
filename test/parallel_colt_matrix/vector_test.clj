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
  (is (= 4 (get-shape v)))
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
  (is (= (pcv/get-vector [2 3]) (subvector v 1 2))))