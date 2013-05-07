;; Testing probabilistic programming scheme with a program that
;; produces a geometric number.

(declare (usual-integrations))

(load "pp")
(load "pp-interface")
(load "constants")
(load "util")
(load "gamma")

(define (accept-reject p* q q-sampler)
  (let* ( (x (q-sampler)) 
          (q-of-x (q x))
          (p*-of-x (p* x))
          (y (uniform-rand 0 q-of-x)))
    (if (< y p*-of-x)
        (begin 
               x)
        (accept-reject p* q q-sampler))))

(define (flip p)
  (< (random 1.0) p))

(define (uniform-rand min max)
  (let ( (x (random (exact->inexact (- max min)))))
    (+ min x)))

(define (uniform-pdf x min max)
  (if (and (<= x max) (>= x min))
      (/ 1 (- max min))
      (error "UNIFORM-PDF: x must be between minimum and maximum")))

(define (beta-pdf x a b)
  (if (and (<= x 1) (>= x 0))
      (* (flo:expt x (- a 1)) (flo:expt (- 1 x) (- b 1)))
      (error "BETA-PDF: x must be between 0 and 1." )))

(define (beta-max a b)
  (define x (/ (- a 1) (- (+ b a) 2)))
  (beta-pdf x a b))

(define (sampler-mc-estimate sampler f niter)
  (define (go n sum)
    (if (> n 0)
        (let ( (x (sampler)))
          (go (- n 1) (+ (f x) sum)))
        (/ sum niter)))
  (go niter 0))

(define (identity x) x)

(define (geometric a b)
  (define p (sample-beta a b))
  (define (go n)
    (if (flip p)
        (go (+ 1 n))
        n))
  (go 0))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; defining pp interface for beta ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; parameters
(define (beta:make-params a b)
  (make-flo-vector (cons a b)))
(define (beta:var params)
  (let ((a (flo:vector-ref params 0))
        (b (flo:vector-ref params 1)))
    (/ (* a b)
       (* (+ a b )
          (+ a b)
          (+ a b 1)))))

;; loglikelihood
(define (beta:log-likelihood val params)
  (let ((a (flo:vector-ref params 0))
        (b (flo:vector-ref params 1))
        (x val))
    (+ (log (gamma (+ a b)))
       (- (log (gamma a)))
       (- (log (gamma b)))
       (* (- a 1) (log x))
       (* (- b 1) (log (- 1 x))))))

;; beta random variates 

(define (beta:rvs params)
  (define a (flo:vector-ref params 0))
  (define b (flo:vector-ref params 1))
  (define (p* x) (beta-pdf x a b))
  (define (q x) (* (uniform-pdf x 0 1) (beta-max a b)))
  (define (q-sampler) (uniform-rand 0 1))
  (accept-reject p* q q-sampler))

  

;; logit
(define (logit x)
  (log (/ x (- 1 x))))

(define (logistic x)
  (/ (exp x) (+ 1 (exp x))))

(define (beta a b #!optional proposer)
  (if (default-object? proposer)
      (set! proposer (proposals:logistic-gaussian (/ (beta:var (beta:make-params  a b))  4))))
  
  (sample 
   'beta
   beta:rvs
   beta:log-likelihood
   (beta:make-params a b)
   proposer))

(define ((proposals:logistic-gaussian var) choice)

  ;; the log likelihood associated with transitioning from x to
  ;; x-prime, where x, x-prime are between 0 and 1. 
  (define (logistic-gaussian:score x x-prime var)
    (let* ( (y       (logit x))
            (y-prime (logit x-prime))
            (ll (gaussian:log-likelihood (- y y-prime)
                                         (gaussian:make-params 0 var)))
            (neg-log-det-jacobian 
             (+ (log x) (log x-prime) (log (- 1 x)) (log (- 1 x-prime)))))
;;      (bkpt '() '())
      (flo:- ll neg-log-det-jacobian)))
    
  (let* ((params (gaussian:make-params 0 var))
         (old-val (choice:val choice))
         (nudge (gaussian:rvs params))
         (new-val (logistic (flo:+ (logit old-val) nudge)))
         (proposal-score (logistic-gaussian:score old-val new-val var)))
    (set! *forward-score* proposal-score)
    (set! *backward-score* proposal-score)
    new-val))
  
;; test
(define (betatest)
  (let ((x (beta 5 1)))
    (emit x 0.2 (likelihood:additive-gaussian 0 0.01))
    x))

(define (geotest)
  (define (go p)
    (display p)(newline)
    (if (flip p)
        (+ 1 (go p))
        0))
  (let* ((p (beta 2 2))
         (count (go p)))
    (emit count 11 likelihood:exact)
    p))
        


  


