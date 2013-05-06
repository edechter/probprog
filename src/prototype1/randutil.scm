(declare (usual-integrations))

;; NOTE: functions in this file use flo-typed math!

(load "constants")
(load "util")

;;;;;;;;;;;;;;
;; Gaussian ;;
;;;;;;;;;;;;;;

;; parameters

(define (gaussian:make-params mean var)
  (vector mean var))
(define (gaussian:mean params)
  (vector-ref params 0))
(define (gaussian:var params)
  (vector-ref params 1))

;; sampling

(define (gaussian:rvs params)
  (let ((mean (gaussian:mean params))
        (var (gaussian:var params)))
    (let ((u (random 1.0))
          (v (random 1.0))
          (std (exact->inexact (sqrt var))))
      (flo:+ mean
             (flo:* std
                    (flo:*
                      (flo:sqrt (flo:* -2. (flo:log u)))
                      (flo:cos (flo:* (flo:* 2. pi) v))))))))

;; log-likelihood

(define (gaussian:log-likelihood val params)
  (let* ((mean (vector-ref params 0))
         (var (vector-ref params 1))
         (centered (exact->inexact (- val mean))))
    (flo:- (flo:/ (flo:* centered centered) (flo:* -2. var))
           (flo:/ (flo:log (flo:* (flo:* 2. pi) var)) 2.))))

;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Discrete/categorical ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;; parameters
(define (discrete:make-params items #!optional weights)
  (let ((tot #f))
    (if (not (default-object? weights))
      (begin
        (set! weights (apply vector weights))
        (set! tot (sigma (lambda (k) (vector-ref weights k)) 0 (fix:- (vector-length weights) 1)))))
    `#(,(list->vector items) ,weights ,tot)))

(define (discrete:items params)
  (vector-ref params 0))

(define (discrete:weights params)
  (vector-ref params 1))

(define (discrete:tot params)
  (vector-ref params 2))

;; sampling

;; if all weights are zero, always chooses first
(define (discrete:rvs params)
  (let ((items (discrete:items params))
        (weights (discrete:weights params))
        (tot (exact->inexact (discrete:tot params))))
    (if (default-object? weights)
      (vector-ref items (random (vector-length items)))
      (let lp ((val (flo:* tot (random 1.0)))
               (idx 0))
        (let ((val (- val (vector-ref weights idx))))
          (if (not (flo:> val 0.))
            (vector-ref items idx)
            (lp val (+ idx 1))))))))

;;NOTE doesn't actually normalize
(define (discrete:log-likelihood val params)
  (let ((items (discrete:items params))
        (weights (discrete:weights params))
        (tot (discrete:tot params)))
    (if (= (vector-length items) 1)
      0.
      (if (default-object? weights)
        (flo:negate (log (vector-length items)))
        (let lp ((idx 0))
          (if (eq? val (vector-ref items idx))
            (log (/ (vector-ref weights idx) tot))
            (lp (+ idx 1))))))))

