(declare (usual-integrations))

;; this file contains the interface functions users can call

(load "randutil")

;;;;;;;;;;;;;;
;; SAMPLING ;;
;;;;;;;;;;;;;;

;; discrete

(define (discrete items #!optional weights proposer)
  (if (and (default-object? proposer)
           (not (default-object? weights))
           (procedure? weights))
    (let ((real-proposer weights))
      (set! weights proposer)
      (set! proposer real-proposer)))

  (sample
    'discrete
    discrete:rvs
    discrete:log-likelihood
    (discrete:make-params items weights)
    proposer))

;; continuous

(define (gaussian mean var #!optional proposer)
  (if (default-object? proposer)
    (set! proposer (proposals:additive-gaussian 0 (/ var 4))))

  (sample
    'gaussian
    gaussian:rvs
    gaussian:log-likelihood
    (gaussian:make-params mean var)
    proposer))

;;;;;;;;;;;;;;
;; EMITTING ;;
;;;;;;;;;;;;;;

(define ((likelihood:additive-gaussian mean var) x obs)
  (gaussian:log-likelihood (exact->inexact obs) (gaussian:make-params (+ mean x) var)))

(define (likelihood:exact x obs)
  (if (eq? x obs) 0. neginf))

;;;;;;;;;;;;;;;
;; PROPOSALS ;;
;;;;;;;;;;;;;;;

(define ((proposals:additive-gaussian mean var) choice)
  (let* ((params (gaussian:make-params mean var))
         (old-val (exact->inexact (choice:val choice)))
         (nudge (gaussian:rvs params))
         (new-val (+ old-val nudge))
         (proposal-score (gaussian:log-likelihood nudge params)))
    (set! *forward-score* proposal-score)
    (set! *backward-score* proposal-score)
    new-val))

(define ((proposals:prior-proposer sampler log-likelihood parameters) choice)
  (let ((new-val (sampler parameters))
        (old-val (choice:val choice)))
    (let ((forward-score (log-likelihood new-val parameters))
          (backward-score (log-likelihood old-val parameters)))
      (set! *forward-score* forward-score)   ;; forward means alternative -> current
      (set! *backward-score* backward-score) ;; backward means current -> alternative
      new-val)))

(define ((proposals:flip p) choice)
  (if (< (random 1.0) p)
      (begin
        (set! *forward-score (flo:log p))
        (set! *backward-score (flo:log (- 1 p)))
        (not (choice:val choice)))
      (begin
        (set! *forward-score (flo:log (- 1 p)))
        (set! *backward-score (flo:log (p)))
        ((choice:val choice)))))
