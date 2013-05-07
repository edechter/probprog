;; A generic API for random variables in a probabilistic program.

(declare (usual-integrations))

(load "pp.scm")

;;;;;;;;;;;;;;;;;;;;;;;;;;
;; rv datatype ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;; value of an rv that has not been sampled
(define the-unsampled-value 'the-unsampled-value)

;; datatype for rv
(define-record-type <rv>
  (rv:new type val params sampler-fn likelihood-fn force-hook force-set-hook proposer)
  rv?
  (type rv:type)
  (val rv:val rv:set-val!)
  (params rv:params rv:set-params!)
  (sampler-fn rv:sampler-fn)
  (likelihood-fn rv:likelihood-fn)
  (force-hook rv:force-hook)
  (force-set-hook rv:force-set-hook)
  (proposer rv:proposer))

(define (rv:type? r type)
  (and (rv? r)
       (eq? (rv:type r) type))) 

(define (rv:unsampled? r)
  (and (rv? r)
       (eq? (rv:val r) the-unsampled-value)))

(define (rv:sampled? r)
  (and (rv? r)
       (not (rv:unsampled? r))))

(define (rv:force r)
  (if (not (rv? r))
      r
      (begin 
        (if (rv:unsampled? r)
            (begin 
              (rv:set-val! r ((rv:sampler-fn r) (rv:params r)))
              ((rv:force-hook r) r)))
        (rv:val r))))

(define (force-set! r val)
  (rv:set-val! r val)
  ((rv:force-set-hook r) r)
  (rv:val r))

(define generic:emit
  (make-generic-procedure (make-procedure-arity 2 #f) 'generic:emit))

(set-generic-procedure-default-generator!
 generic:emit (lambda (proc tags) emit))

