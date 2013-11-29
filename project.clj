(defproject parallel-colt-matrix "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [net.mikera/core.matrix "0.16.0"]
                 [net.sourceforge.parallelcolt/parallelcolt "0.10.0"]
                 [org.clojure/tools.nrepl "0.2.0"]]
  :profiles {:dev                                                               
             {:dependencies                                                     
              [[criterium "0.3.1"] ;; Benchmarking
               ]}}
  :warn-on-reflection true)
