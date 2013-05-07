

(declare (usual-integrations))

(load "pp")
(load "pp-interface")
(load "pp-records")

(load "generic-rv")
(load "beta")
(load "bern")

;; a global beta-bernoulli data structure that keeps track of all
;; bernoullis conjugate to the beta
(define-record-type <beta-bern-ds>
  (beta-bern-ds:new a b berns)
  beta-bern-ds?
  (a beta-bern-ds:a beta-bern-ds:set-a!)
  (b beta-bern-ds:b beta-bern-ds:set-b!)
  (berns beta-bern-ds:berns beta-bern-ds:set-berns!))

(define (beta-bern-ds:init)
  (beta-bern-ds:new 'no-a 'no-b (make-eq-hash-table)))

(define (beta-bern-ds:insert-bern berns bern val)
  (hash-table/put! berns bern val))

(define *beta-bern-ds* (beta-bern-ds:init))

(define (*beta-bern-ds*:show)
  (display "a: ") (display (beta-bern-ds:a *beta-bern-ds*)) (newline)
  (display "b: ") (display (beta-bern-ds:b *beta-bern-ds*)) (newline)
  (display "berns: ") (display (hash-table->alist (beta-bern-ds:berns *beta-bern-ds*))) (newline))
  

(define (*beta-bern-ds*:insert-bern! bern val)
  (beta-bern-ds:insert-bern (beta-bern-ds:berns *beta-bern-ds*)
                            bern 
                            val))

(define (*beta-bern-ds*:update-counts!)
  (let* ( (counts (hash-table/datum-list (beta-bern-ds:berns *beta-bern-ds*)))
          (heads (filter (lambda (x) (eq? 1 x)) counts))
          (tails (filter (lambda (x) (eq? 0 x)) counts)))
    (beta-bern-ds:set-a! *beta-bern-ds* (+ (length heads) (beta-bern-ds:a *beta-bern-ds*)))
    (beta-bern-ds:set-b! *beta-bern-ds* (+ (length tails) (beta-bern-ds:b *beta-bern-ds*)))))
          

(define (beta a b #!optional proposer)
  (if (default-object? proposer)
      (set! proposer (proposals:logistic-gaussian 
                      (/ (beta:var (beta:make-params  a b))  4))))

  (beta-bern-ds:set-a! *beta-bern-ds* a)
  (beta-bern-ds:set-b! *beta-bern-ds* b)

  (let ( (rv (rv:new 'beta
                     the-unsampled-value
                     (pair->flo-vector (cons a b))
                     beta:rvs-global
                     beta:log-likelihood
                     beta:force-hook
                     beta:force-set-hook
                     proposer)))
    rv))

(define (beta:rvs-global params)
  (let ((a (beta-bern-ds:a *beta-bern-ds*))
        (b (beta-bern-ds:b *beta-bern-ds*)))
    (beta:rvs (pair->flo-vector (cons a b)))))

(define (beta:force-hook rv)
  'unspecific)

(define (beta:force-set-hook rv)
  'unspecific)

(define (bern p #!optional proposer)

  (if (default-object? proposer)
      (set! proposer (proposals:flip 0.5)))
  
  (cond ( (number? p)
          (let ((rv (rv:new 'bern
                            the-unsampled-value
                            p
                            bern:rvs 
                            bern:log-likelihood 
                            'unspecific
                            'unspecific
                            proposer)))
            rv))
          ;; if p is an unsampled beta rv
        ((and (rv:type? p 'beta) (rv:unsampled? p))
         (let ((rv (rv:new 'bern
                           the-unsampled-value
                           ;; params
                           p
                           ;; sampler
                           (lambda (beta-rv)
                               (bern:rvs (rv:force beta-rv)))
                           ;; marginal log-likelihood TODO: currently
                           ;; assumes there is only one beta in the
                           ;; world
                           (lambda (x beta-rv)
                             (bern:marginal-log-likelihood x
                                                           (beta-bern-ds:a *beta-bern-ds*)
                                                           (beta-bern-ds:b *beta-bern-ds*)))
                           ;; force-hook
                           bern:force-hook
                           bern:force-set-hook
                           proposer)))
             rv))
          ((and (rv:type? p 'beta) (rv:sampled? p))
           (let ((rv (rv:new 'bern
                           the-unsampled-value
                           ;; params
                           p
                           ;; sampler
                           (lambda (beta-rv)
                               (bern:rvs (rv:val beta-rv)))
                           ;; conditional log-likelihood 
                           (lambda (x beta-rv)
                             (bern:rvs x (rv:val beta-rv)))
                           ;; force-hook
                           bern:force-hook
                           ;; force-set-hook
                           bern:force-set-hook
                           proposer)))
             rv))
          ( (error "BERN: don't know what to do with this type of parameter." p))))

(define (bern:force-hook bern)
  (*beta-bern-ds*:insert-bern! bern (rv:val bern))
  (*beta-bern-ds*:update-counts!))
  
(define (bern:force-set-hook bern )
  (display "Bern force set hook!")
  (*beta-bern-ds*:insert-bern! bern (rv:val bern))
  (*beta-bern-ds*:update-counts!))

(define (emit-bern bern val)
  (force-set! bern val))
  
;; Goal program
(define (test-beta-bern)
  (let* ((p (beta 1 2))
        (x1 (bern p))
        (x2 (bern p))
        (x3 (bern p))
        (x4 (bern p))
        (x5 (bern p))
        (x6 (bern p)))
    (emit-bern x1 0)
    (emit-bern x2 0)
    (emit-bern x3 0)
    (emit-bern x4 0)
    (emit-bern x5 0)
    (emit-bern x6 0)
    p))





  


