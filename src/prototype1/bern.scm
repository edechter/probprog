
(load "gamma")

(define (bern:rvs p) ;; params is singleton vector of a p in [0, 1]
    (< (random 1.0) p))

(define (bern:log-likelihood val p)
  (if val
      (flo:log p) 
      (flo:log (- 1 p))))

  
(define (bern:marginal-log-likelihood x a b)
  (let ((y (if x 1 0)))
    (+ ( log (+ a b))
       ( log (gamma (+ y a)))
       ( log (gamma (+ b (- x) 1)))
       ( - ( log (gamma a)))
       ( - ( log (gamma b))))))
